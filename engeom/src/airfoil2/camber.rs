use crate::AngleDir::Ccw;
use crate::airfoil2::SectionInput;
use crate::airfoil2::inscribed::{Inscribed, InscribedVec};
use crate::common::dist;
use crate::geom2::hull::{convex_hull_2d, farthest_pair_on_hull};
use crate::geom2::{LineOps2, rot90};
use crate::{Curve2, Line2, Result, SurfacePoint2};

const EDGE_STOP_FRACTION: f64 = 0.375;

/// Extract the unambiguous inscribed circles of an airfoil section.
///
/// - The circles will be returned in a consistent order, but the order may be _either_ from leading
///   to trailing edge _or_ the reverse.
/// - The `p0` and `p1` points will all be oriented in a consistent way, but the order may be
///   _either_ from upper to lower surface _or_ the reverse.
///
/// The camber line extraction algorithm will terminate when the distance left to the farthest
/// point on the edge is more than 3/8's of the radius of the last inscribed circle beyond its
/// perimeter.  From there, edge-aware algorithms should be used to extend the camber line and any
/// remaining inscribed circles.
///
/// # Arguments
///
/// * `section`: An airfoil section represented by a curve. May be open at the leading or trailing
///   edge but not both.  Should not be open in the middle unless the gap is small compared to the
///   thickness of the airfoil at the open portion.
/// * `tol`: A tolerance for extracting the circles. More circles will be added until interpolation
///   between neighboring circles has an error less than this amount. The circle centers will be
///   located with a tolerance 1/10th of this value.
///
/// returns: Result<Vec<Inscribed, Global>, Box<dyn Error, Global>>
pub fn extract_inscribed_circles(input: &SectionInput) -> Result<Vec<Inscribed>> {
    let hull = convex_hull_2d(input.section.points());
    let (i0, i1) = farthest_pair_on_hull(input.section.points(), &hull);
    let naive_chord = Line2::from_points(&input.section.points()[i0], &input.section.points()[i1]);

    // The starting line will be halfway down, but ultimately we're going to want to be able to
    // sample at multiple places along the chord and look for valid crossing lines, with some
    // mechanism to guess at their quality
    //
    // Alternately, perhaps something with the zhang suen line can be a good starting point
    let ref_line = Line2::new(naive_chord.at(0.5), naive_chord.orthogonal());

    // First we'll gather the circles in the front half
    let start_line = input.crossing_line(&ref_line).ok_or_else(|| {
        "Failed to find a crossing line for the initial reference line".to_string()
    })?;
    let mut working = extract_half_circles(input, &start_line)?;

    // We'll reverse the vector so that we can add the elements from the back half of the section,
    // removing the last element since it will be the same as the first element of the back half.
    working.reverse_order();
    working.pop();

    // Now we gather the circles in the back half
    let start_line = input.crossing_line(&ref_line.reversed()).ok_or_else(|| {
        "Failed to find a crossing line for the reversed reference line".to_string()
    })?;
    working.extend(extract_half_circles(input, &start_line)?);

    Ok(working.take_vec())
}

/// Extracts the inscribed circles starting at the initial line and proceeding in the direction of
/// the line orthogonal direction until reaching the edge neighborhood or throwing an error.
///
/// # Arguments
///
/// * `section`: An airfoil section represented by a `Curve2`
/// * `line`: A valid crossing line for the airfoil (t=0 and t=1 must be on the section perimeter)
/// * `tol`: A tolerance for extracting the circles. More circles will be added until interpolation
///   between neighboring circles has an error less than this amount. The circle centers will be
///   located with a tolerance 1/10th of this value.
///
/// returns: Result<Vec<Inscribed, Global>, Box<dyn Error, Global>>
fn extract_half_circles(input: &SectionInput, line: &Line2) -> Result<InscribedVec> {
    let mut results = InscribedVec::empty();
    let mut working_line = *line;

    loop {
        let circle = input.inscribed_from_crossing(&working_line);
        results.refine_and_push(circle, input);
        let last = results
            .last()
            .ok_or("Failed to get last inscribed circle".to_string())?;

        match advance_search(input, last)? {
            Some(next_line) => {
                working_line = next_line;
            }
            None => {
                break;
            }
        };
    }

    Ok(results)
}

/// Attempts to advance the inscribed circle search by jumping forward from the last inscribed
/// circle by a fraction of its radius and seeking a valid crossing line. If it fails it will
/// attempt again at reduced distances. If it is close to the edge it will return `None` to indicate
/// that we are getting close to the ambiguous portion of the MCL and the search should terminate.
///
/// # Arguments
///
/// * `section`: An airfoil section represented by a `Curve2`
/// * `last`: The last valid inscribed circle
/// * `tol`: The circle finding tolerance, which is used as a lower limit of when to terminate the
///   advancement attempt when all else fails.
///
/// returns: Result<Option<Line2>, Box<dyn Error, Global>>
///
/// # Warnings
///
/// If the airfoil is curved enough that part of the other side lies in front of the current
/// station's search direction, this function will produce unknown results and will need to be
/// modified.
fn advance_search(input: &SectionInput, last: &Inscribed) -> Result<Option<Line2>> {
    // We will begin by finding the camber point/direction of the last station, which will be used
    // to jump forward and create a new spanning ray.  However, we'll first check the distance from
    // the camber point to the farthest point on the section in the camber direction.  As we get
    // closer to the edge of the airfoil, we will want to terminate the search.
    let cp = last.camber_point();

    // We unwrap this because the only way it would fail is if the section is empty, which
    // would have prevented us from getting here in the first place.
    // let (_, farthest) = section.max_point_in_direction(&cp.normal).unwrap();
    // let distance = cp.scalar_projection(&farthest);
    let distance = remaining(input.section, &cp);

    // When the distance beyond the last inscribed circle is less than 3/8 of the circle's radius,
    // we will consider ourselves close enough to the edge of the airfoil to terminate the search.
    // Getting closer to the edge will increase the probability that the assumptions of no local
    // maxima along the ray are violated.
    if distance - last.c.r() < last.c.r() * EDGE_STOP_FRACTION {
        return Ok(None);
    }

    // Now we will create a new spanning ray which will be used to find the next inscribed circle.
    // We will start by jumping forward 1/8 of the last circle's radius, and we will adjust this
    // value down as we have failures.  So long as we move forward at least 5% of the last circle's
    // radius, we will consider the search to have advanced.
    let mut frac = 0.125;
    while frac > 0.05 {
        let advance_dist = frac * last.c.r();
        if advance_dist < input.general_tol {
            return Ok(None);
        }
        let next_center = cp.at_distance(advance_dist);
        let test_dir = rot90(Ccw) * cp.normal;
        let test_line = Line2::new(next_center, test_dir.into_inner());

        if let Some(line) = input.crossing_line(&test_line) {
            // First, we want to test if the new ray spans at least 50% of the last station's
            // distance between the upper and lower contact points.  This is a heuristic to ensure
            // we haven't taken a step where the section thickness is dropping off too quickly.
            let last_dist = dist(&last.p0, &last.p1);
            let new_dist = line.direction.norm();

            if new_dist < 0.5 * last_dist {
                frac *= 0.75;
                continue;
            }

            return Ok(Some(line));
        } else {
            frac *= 0.75;
        }
    }

    Err("Failed to advance search".into())
}

/// This is a local helper function which tries to find the maximum remaining distance in front of
/// a surface point. It's used to see how close to the end of the section we are.
fn remaining(section: &Curve2, cp: &SurfacePoint2) -> f64 {
    let mut result = f64::NEG_INFINITY;
    for p in section.points() {
        if cp.scalar_projection(p).is_sign_positive() {
            result = result.max(dist(cp, p));
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve2;
    use crate::tests::airfoil_curve;
    use approx::assert_abs_diff_eq;

    fn len_dist(curve: &Curve2, inscribed: &Inscribed) -> (f64, f64) {
        let m = curve.at_closest_to_point(&inscribed.c);
        let distance = dist(&m, &inscribed.c);
        let len = m.length_along();
        (len, distance)
    }

    #[test]
    fn simple_extract() -> Result<()> {
        let section = airfoil_curve();
        let input = SectionInput::new(&section, 1e-3);
        let mut circles = extract_inscribed_circles(&input)?;

        #[rustfmt::skip]
        let expected = curve2!((-10.9929809799, -10.0096070029), (-11.1856809788, -8.4210146505), (-11.1979663503, -8.3009068172), (-11.2073856723, -8.1805064424), (-11.2300944957, -7.9401148484), (-11.2370146936, -7.8801337205), (-11.2422757600, -7.8199838272), (-11.2533002501, -7.6997293900), (-11.2605164039, -7.5792258431), (-11.2641093798, -7.4585108797), (-11.2705592482, -7.2170527616), (-11.2723924681, -7.1567008013), (-11.2713028815, -7.0963316450), (-11.2690011327, -6.9755925932), (-11.2623769903, -6.7342371419), (-11.2516253651, -6.4930243979), (-11.2407717504, -6.3171234209), (-11.2272524867, -6.1414159317), (-11.1948367499, -5.7695101241), (-11.1744259159, -5.5733254593), (-11.1503741160, -5.3775228415), (-11.1243875016, -5.1712653215), (-11.0943079427, -4.9655693012), (-11.0626053811, -4.7489535259), (-11.0468077204, -4.6406370847), (-11.0286172745, -4.5327024053), (-10.9467566114, -4.0800271140), (-10.9014017434, -3.8430237307), (-10.8537543864, -3.6064661566), (-10.8276799476, -3.4825520492), (-10.7995119907, -3.3591007079), (-10.7417682385, -3.1125384615), (-10.6761724981, -2.8551276971), (-10.6074481243, -2.5985243962), (-10.5322841724, -2.3302585809), (-10.4510406333, -2.0637324727), (-10.3623155273, -1.7854191911), (-10.2656144644, -1.5097655429), (-10.1602831879, -1.2224815515), (-10.1050979014, -1.0797957040), (-10.0474140414, -0.9380898440), (-9.9862900239, -0.7900654609), (-9.9232832989, -0.6428289765), (-9.8593791703, -0.4959876662), (-9.7933700287, -0.3500829824), (-9.7214698048, -0.1992337287), (-9.6467081582, -0.0497613842), (-9.5722986519, 0.0998898048), (-9.4959688732, 0.2485612745), (-9.4558483273, 0.3257167612), (-9.4139625062, 0.4019242919), (-9.3293375190, 0.5538740917), (-9.1514116343, 0.8526782192), (-9.0526360338, 1.0029533154), (-8.9510771512, 1.1513804244), (-8.8469239049, 1.2979614361), (-8.7365578367, 1.4398907691), (-8.6199032273, 1.5822686653), (-8.4999477838, 1.7218719604), (-8.4389218620, 1.7907205731), (-8.3762627655, 1.8580583860), (-8.2457209651, 1.9877139624), (-8.1785271080, 2.0516114248), (-8.1089536789, 2.1129065019), (-8.0379038247, 2.1724218456), (-7.9643922128, 2.2288567040), (-7.8887466289, 2.2824817244), (-7.8118376707, 2.3343077113), (-7.7338414455, 2.3844993076), (-7.6538982920, 2.4315932052), (-7.5748766242, 2.4754833402), (-7.4939490505, 2.5157468938), (-7.4118102841, 2.5536563592), (-7.3285325676, 2.5889917723), (-7.2443629547, 2.6222740282), (-7.1594502149, 2.6536708702), (-6.9892963392, 2.7154521102), (-6.6559944164, 2.8196774943), (-6.4881178402, 2.8674769092), (-6.3193578084, 2.9119323324), (-5.9906905832, 2.9863852591), (-5.6603553579, 3.0529203721), (-5.3404088010, 3.1098973737), (-5.1800093938, 3.1358479366), (-5.0192377385, 3.1594793121), (-4.8641801280, 3.1820960343), (-4.7088364682, 3.2025924415), (-4.3977112667, 3.2402194826), (-4.0971828004, 3.2739836627), (-3.9467344249, 3.2890761048), (-3.7960902567, 3.3021080193), (-3.2143623637, 3.3457186102), (-2.9332627519, 3.3609284880), (-2.6520226878, 3.3736803942), (-2.3803733679, 3.3813929657), (-2.1086562587, 3.3862365938), (-1.8462608130, 3.3867070478), (-1.5838783955, 3.3848825045), (-1.0770516568, 3.3724511155), (-0.8320468104, 3.3625058674), (-0.5871438920, 3.3501701750), (-0.3501791799, 3.3368901824), (-0.1133416886, 3.3215489894), (0.1159430298, 3.3055089204), (0.3450645073, 3.2873007743), (0.5671589583, 3.2691815710), (0.7890562313, 3.2488609393), (1.2191426586, 3.2049076034), (1.6360889536, 3.1560672521), (2.0406192053, 3.1046645890), (2.4330734035, 3.0498901687), (2.8142818184, 2.9943324767), (3.1841337563, 2.9336024832), (3.5534833957, 2.8699042457), (3.9126733060, 2.8071827759), (4.2612455166, 2.7412738103), (4.6003651054, 2.6762603466), (4.9302204296, 2.6108104939), (5.2510679756, 2.5453756711), (5.5630763674, 2.4788372525), (5.8670531599, 2.4137140639), (6.1628857951, 2.3486071114), (6.4506201959, 2.2817737249), (6.7308874555, 2.2146542088), (7.0043834612, 2.1484523194), (7.2709510704, 2.0816262449), (7.5312571512, 2.0155029005), (7.7853046546, 1.9493797564), (8.0334163425, 1.8834960893), (8.2757983586, 1.8176550542), (8.5129960695, 1.7530692590), (8.7449340818, 1.6891456137), (8.9715777688, 1.6249749899), (9.1933507775, 1.5615887998), (9.4101977297, 1.4983037156), (9.6224967028, 1.4355515943), (9.8303738159, 1.3733913992), (10.0339670692, 1.3119948429), (10.2334392353, 1.2514015282), (10.4286763930, 1.1909259471), (10.6201151025, 1.1314476533), (10.8078335266, 1.0728663779), (10.9916521529, 1.0144933724), (11.1717652802, 0.9566685723), (11.3484462702, 0.8997658135), (11.5217574219, 0.8438318016), (11.6916207185, 0.7883938762), (11.8582446193, 0.7337162600), (12.0217540779, 0.6798890872), (12.1819937496, 0.6263063298), (12.3392720652, 0.5734982494), (12.4936560083, 0.5213797273), (12.6451271584, 0.4696408183), (12.7938768151, 0.4185944835), (12.9399536758, 0.3681781361), (13.0834883348, 0.3184971764), (13.2244088981, 0.2691421281), (13.3630071589, 0.2206287232), (13.4993571678, 0.1730168599), (13.6333092147, 0.1257194486), (13.7650410757, 0.0790061288), (13.8946163307, 0.0328655851), (14.0220843746, -0.0126431598), (14.1474647870, -0.0575821357), (14.2709206553, -0.1017786143), (14.3925501705, -0.1450635516), (14.5122819983, -0.1876314780), (14.6299310815, -0.2300874959), (14.7457813118, -0.2719130770), (14.8598788712, -0.3131610180), (14.9721530929, -0.3542252914), (15.0828288723, -0.3947394378), (15.1919547434, -0.4346464314), (15.2993899147, -0.4744092608), (15.4053239285, -0.5137206378), (15.5097677730, -0.5526802583), (15.6126822754, -0.5915022881), (15.7142842409, -0.6298140873), (15.8145494621, -0.6677630919), (15.9135918457, -0.7052130657), (16.0112620683, -0.7426138745), (16.1077694839, -0.7795704860), (16.2030831975, -0.8162362090), (16.2972583815, -0.8525256381), (16.3902841407, -0.8885628944), (16.4822396770, -0.9242030838), (16.5731970038, -0.9592984999), (16.6631134812, -0.9940125052), (16.7519798520, -1.0284364122), (16.8398664238, -1.0624504253), (16.9267694657, -1.0960729718), (17.0126679172, -1.1293536829), (17.0976375882, -1.1621287627), (17.1816085545, -1.1946120514), (17.2646708670, -1.2266081299), (17.3468087873, -1.2581543170), (17.4280368862, -1.2892007179), (17.5082559970, -1.3200396191), (17.5875667772, -1.3504416199), (17.6658872697, -1.3806373972), (17.7431447997, -1.4108342366), (17.8195492522, -1.4406309049), (17.8950961816, -1.4700656286), (17.9697783581, -1.4991578882), (18.0436494740, -1.5277725316), (18.1166095474, -1.5561966672), (18.1886867078, -1.5843597872), (18.2599047949, -1.6121342582), (18.3302814787, -1.6394877195), (18.3997442260, -1.6666165682), (18.4684002152, -1.6932773696), (18.5361567214, -1.7197396101), (18.6031009390, -1.7458103888), (18.6692373216, -1.7715096854), (18.7345231814, -1.7969777116), (18.7990484207, -1.8220065210), (18.8626368163, -1.8470545549), (18.9253779477, -1.8718046944), (18.9873041213, -1.8961982421), (19.0484317668, -1.9202295080), (19.1086837901, -1.9441180917), (19.1681578911, -1.9676770273), (19.2268498836, -1.9909651460), (19.2848262034, -2.0138569160), (19.3420250286, -2.0365462860), (19.3984799075, -2.0589772036), (19.4542014947, -2.0811980742), (19.5092010046, -2.1031952784), (19.5634011977, -2.1251332445), (19.6169200626, -2.1466921485), (19.6696895511, -2.1680677815), (19.7217605948, -2.1891610312), (19.7730707661, -2.2101314987), (19.8237713567, -2.2306555565), (19.8738359343, -2.2507640768), (19.9230983394, -2.2707619261), (19.9716463603, -2.2905006890), (20.0197002785, -2.3094416785), (20.0669596178, -2.3281642159), (20.1134670406, -2.3466077073), (20.1595742904, -2.3639466160), (20.2052398909, -2.3798931777), (20.2507940439, -2.3932930066), (20.2951863859, -2.4056878341), (20.7304043107, -2.5021209041))?;

        // Do a quick orientation just so that we know we're increasing
        let (l0, _) = len_dist(&expected, &circles[0]);
        let (l1, _) = len_dist(&expected, &circles[circles.len() - 1]);
        if l1 < l0 {
            circles.reverse();
        }

        let mut last_len = f64::NEG_INFINITY;
        for c in circles.iter() {
            let (l, d) = len_dist(&expected, c);
            assert!(
                l > last_len,
                "Circle positions along the expected mcl are not increasing"
            );
            assert_abs_diff_eq!(d, 0.0, epsilon = 1e-3);
            last_len = l;
        }

        Ok(())
    }
}
