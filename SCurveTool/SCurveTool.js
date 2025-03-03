const M_PI = Math.PI;
const M_2PI = M_PI * 2.0;
const GRAVITY_MSS = 9.81;

pos_plot = {}
vel_plot = {}
accel_plot = {}
jerk_plot = {}
wp_pos_plot = {}
function initial_load()
{
    const time_scale_label = "Time (s)";
    let plot;

    // Waypoints
    wp_pos_plot.data = [{ type:'scatter3d',  x:[], y:[], z:[], name: 'WP', mode: 'lines+markers', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s³" },
                        { type:'scatter3d',  x:[], y:[], z:[], name: 'Target', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s³", line: { width: 10 }}];

    wp_pos_plot.layout = {
        legend: { itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 60, r: 50, t: 20 },
        scene: {
            xaxis: { title: {text: "East (m)" }, autorange: "reversed" },
            yaxis: { title: {text: "North (m)"}, autorange: "reversed" },
            zaxis: { title: {text: "Up (m)" } }
        }
    }

    plot = document.getElementById("waypoint_plot")
    Plotly.purge(plot)
    Plotly.newPlot(plot, wp_pos_plot.data, wp_pos_plot.layout, { displaylogo: false })

    // Jerk
    jerk_plot.data = [{ x:[], y:[], name: 'scurve-log', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s³" },
                      { x:[], y:[], name: '---', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s³" }];

    jerk_plot.layout = {
        legend: { itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 60, r: 50, t: 20 },
        xaxis: { title: {text: time_scale_label } },
        yaxis: { title: {text: "Jerk (m/s³)" } },
        showlegend: true
    }

    plot = document.getElementById("jerk_plot")
    Plotly.purge(plot)
    Plotly.newPlot(plot, jerk_plot.data, jerk_plot.layout, { displaylogo: false })

    // Acceleration
    accel_plot.data = [{ x:[], y:[], name: 'scurve-log', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s²" },
                       { x:[], y:[], name: 'output', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s²" }]

    accel_plot.layout = {
        legend: { itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 60, r: 50, t: 20 },
        xaxis: { title: {text: time_scale_label } },
        yaxis: { title: {text: "Acceleration (m/s²)" } },
        showlegend: true
    }

    plot = document.getElementById("accel_plot");
    Plotly.purge(plot);
    Plotly.newPlot(plot, accel_plot.data, accel_plot.layout, { displaylogo: false });

    // velocity
    vel_plot.data = [{ x:[], y:[], name: 'scurve-log', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s" },
                     { x:[], y:[], name: 'output', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s" }];

    vel_plot.layout = {
        legend: { itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 60, r: 50, t: 20 },
        xaxis: { title: {text: time_scale_label } },
        yaxis: { title: {text: "Velocity (m/s)" } },
        showlegend: true,
        shapes: [{
            type: 'line',
            line: { dash: "dot" },
            xref: 'paper',
            x0: 0,
            x1: 1,
            visible: false,
        }]
    }

    plot = document.getElementById("vel_plot")
    Plotly.purge(plot)
    Plotly.newPlot(plot, vel_plot.data, vel_plot.layout, { displaylogo: false })

    // position
    pos_plot.data = [{ x:[], y:[], name: 'scurve-log', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m" },
                     { x:[], y:[], name: 'output', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m" }]
                    //  { x:[], y:[], name: 'Y', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m" },
                    //  { x:[], y:[], name: 'Z', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m" }]

    pos_plot.layout = {
        legend: { itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 60, r: 50, t: 20 },
        xaxis: { title: {text: time_scale_label } },
        yaxis: { title: {text: "Position (m)" } },
        showlegend: true,
        shapes: [{
            type: 'line',
            line: { dash: "dot" },
            xref: 'paper',
            x0: 0,
            x1: 1,
            visible: false,
        }]
    }

    plot = document.getElementById("pos_plot")
    Plotly.purge(plot)
    Plotly.newPlot(plot, pos_plot.data, pos_plot.layout, { displaylogo: false })


    // Link all time axis
    link_plot_axis_range([
        ["jerk_plot", "x", "", jerk_plot],
        ["accel_plot", "x", "", accel_plot],
        ["vel_plot", "x", "", vel_plot],
        ["pos_plot", "x", "", pos_plot],
    ])

    // // Link plot reset
    link_plot_reset([
        ["jerk_plot", jerk_plot],
        ["accel_plot", accel_plot],
        ["vel_plot", vel_plot],
        ["pos_plot", pos_plot],
    ])
}

// Utility functions
const FLT_EPSILON = 1.1920928955078125e-7;
function is_zero(value)
{
    return (Math.abs(value) < FLT_EPSILON);
}

function is_positive(value)
{
    return value >= FLT_EPSILON;
}

function is_negative(value)
{
    return value <= -1.0 * FLT_EPSILON;
}

function is_equal(v1, v2)
{
    return is_zero(v1 - v2);
}

function radians(deg)
{
    return deg * M_PI / 180.0;
}

// convert a maximum lean angle in degrees to an accel limit in m/s/s
function angle_to_accel(angle_deg)
{
    return GRAVITY_MSS * Math.tan(radians(angle_deg));
}

// kinematic_limit calculates the maximum acceleration or velocity in a given direction.
// based on horizontal and vertical limits.
function kinematic_limit(direction, max_xy, max_z_pos, max_z_neg)
{
    if (is_zero(direction.length_squared()) || is_zero(max_xy) || is_zero(max_z_pos) || is_zero(max_z_neg)) {
        return 0.0;
    }

    max_xy = Math.abs(max_xy);
    max_z_pos = Math.abs(max_z_pos);
    max_z_neg = Math.abs(max_z_neg);

    let unit_direction = direction.normalize();
    const xy_length = new Vector(unit_direction.x, unit_direction.y, 0.0).length();

    if (is_zero(xy_length)) {
        // only moving in z
        return is_positive(unit_direction.z) ? max_z_pos : max_z_neg;
    }

    if (is_zero(unit_direction.z)) {
        // only moving in xy
        return max_xy;
    }

    const slope = unit_direction.z/xy_length;
    if (is_positive(slope)) {
        if (Math.abs(slope) < max_z_pos/max_xy) {
            return max_xy/xy_length;
        }
        return Math.abs(max_z_pos/unit_direction.z);
    }

    // Got this far then we are looking at a descending slope
    if (Math.abs(slope) < max_z_neg/max_xy) {
        return max_xy/xy_length;
    }
    return Math.abs(max_z_neg/unit_direction.z);
}


const SegmentType = Object.freeze({
    CONSTANT_JERK: 0,
    POSITIVE_JERK: 1,
    NEGATIVE_JERK: 2
});

class Segments {
    constructor()
    {
        this.jerk_ref = 0;     // jerk reference value for time segment (the jerk at the beginning, middle or end depending upon the segment type)
        this.seg_type = 0;     // segment type (jerk is constant, increasing or decreasing)
        this.end_time = 0;     // final time value for segment
        this.end_accel = 0;    // final acceleration value for segment
        this.end_vel = 0;      // final velocity value for segment
        this.end_pos = 0       // final position value for segment
    }
}

class KinematicLogging {
    constructor()
    {
        this.time = [];     // (s) Array of time
        this.pos = [];      // (m) Array of type Vector
        this.vel = [];      // (m/s) Array of type Vector
        this.accel = [];    // (m/s/s) Array of type Vector
        this.jerk = [];     // (m/s/s/s) Array of type Vector
    }
}

class SCurve {
    constructor() {

        this.segments_max = 23; // maximum number of time segments

        // init the array of empty segments
        this.segment = Array.from({ length: this.segments_max }, (_, i) => new Segments());

        this.track = new Vector();
        this.delta_unit = new Vector();

        this.SEG_INIT = 0;
        this.SEG_ACCEL_END = 7;
        this.SEG_SPEED_CHANGE_END = 14;
        this.SEG_CONST = 15;
        this.SEG_TURN_OUT = 15;
        this.SEG_TURN_IN = 4;
        this.SEG_DECEL_END = 22;

        this.logger = new KinematicLogging();
    }

    // Initialization function
    init() {
        this.snap_max = 0.0;
        this.jerk_max = 0.0;
        this.accel_max = 0.0;
        this.vel_max = 0.0;
        this.time = 0.0;  // time that defines position on the path

        this.num_segs = 0;
        this.num_segs = this.add_segment(this.num_segs, 0.0, SegmentType.CONSTANT_JERK, 0.0, 0.0, 0.0, 0.0)

        this.track.zero();
        this.delta_unit.zero();
        this.position_sq = 0.0;
    }

    // increment the internal time
    advance_time(dt) {
        this.time = Math.min(this.time+dt, this.time_end());
    }

    // time at the end of the sequence
    time_end()
    {
        if (this.num_segs != this.segments_max) {
            return 0.0;
        }
        return this.segment[this.SEG_DECEL_END].end_time;
    }

    // Calculate track motion profile between two 3D points
    calculate_track(origin, destination, speed_xy, speed_up, speed_down, accel_xy, accel_z, snap_maximum, jerk_maximum) {
        this.init();

        // Compute vector between origin and destination
        let track_temp = destination.subtract(origin);
        if (track_temp.is_zero() || is_zero(track_temp.length_squared())) {
            return;
        }

        this.snap_max = snap_maximum;
        this.jerk_max = jerk_maximum;

        this.set_kinematic_limits(origin, destination, speed_xy, speed_up, speed_down, accel_xy, accel_z);

        if (!is_positive(this.snap_max) || !is_positive(this.jerk_max) || 
            !is_positive(this.accel_max) || !is_positive(this.vel_max)) {
            throw new Error("SCurve: Invalid kinematic parameters");
        }

        this.track = track_temp;
        const track_length = this.track.length();

        if (is_zero(track_length)) {
            this.delta_unit.zero();
        } else {
            // with non-zero track length we can build the full s-curve?
            this.delta_unit = this.track.normalize();
            this.add_segments(track_length);
        }

        // Check if all of the segments of the s-curve are valid
        if (!this.is_valid()) {
            throw new Error("SCurve::calculate_track invalid path")
        }
    }

    // Set kinematic limits
    set_kinematic_limits(origin, destination, speed_xy, speed_up, speed_down, accel_xy, accel_z) {
        // ensure arguments are positive
        speed_xy = Math.abs(speed_xy);
        speed_up = Math.abs(speed_up);
        speed_down = Math.abs(speed_down);
        accel_xy = Math.abs(accel_xy);
        accel_z = Math.abs(accel_z);

        let direction = destination.subtract(origin);
        this.vel_max = kinematic_limit(direction, speed_xy, speed_up, speed_down);
        this.accel_max = kinematic_limit(direction, accel_xy, accel_z, accel_z);
    }

    // add single S-Curve segment
    // populate the information for the segment specified in the path by the index variable.
    // the index variable is incremented to reference the next segment in the array
    add_segment(index, end_time, seg_type, jerk_ref, end_accel, end_vel, end_pos)
    {
        this.segment[index].end_time = end_time;
        this.segment[index].seg_type = seg_type;
        this.segment[index].jerk_ref = jerk_ref;
        this.segment[index].end_accel = end_accel;
        this.segment[index].end_vel = end_vel;
        this.segment[index].end_pos = end_pos;
        return index + 1;
    }


    // generate the segments for a path of length L
    // the path consists of 23 segments
    // 1 initial segment
    // 7 segments forming the acceleration S-Curve
    // 7 segments forming the velocity change S-Curve
    // 1 constant velocity S-Curve
    // 7 segments forming the deceleration S-Curve
    add_segments(L)
    {
        if (is_zero(L)) {
            return;
        }

        let [Jm, tj, t2, t4, t6] = this.calculate_path(this.snap_max, this.jerk_max, 0.0, this.accel_max, this.vel_max, L * 0.5);

        // Acceleration s-curve phase
        this.num_segs = this.add_segments_jerk(this.num_segs, tj, Jm, t2);    // +3 segs (4) - increasing acceleration phase (positive jerk)
        this.num_segs = this.add_segment_const_jerk(this.num_segs, t4, 0.0);  // +1 seg (5) - constant acceleration (zero jerk)
        this.num_segs = this.add_segments_jerk(this.num_segs, tj, -Jm, t6);   // +3 segs (8) - decreasing acceleration phase (negative jerk)

        // remove numerical errors
        // We expect that we will have reached the 0 acceleration, constant velocity phase by this point
        this.segment[this.SEG_ACCEL_END].end_accel = 0.0;

        // add empty speed adjust segments
        this.num_segs = this.add_segment_const_jerk(this.num_segs, 0.0, 0.0); // +1 seg - velocity ??
        this.num_segs = this.add_segment_const_jerk(this.num_segs, 0.0, 0.0); // +1 seg - velocity ??
        this.num_segs = this.add_segment_const_jerk(this.num_segs, 0.0, 0.0); // +1 seg - velocity ??
        this.num_segs = this.add_segment_const_jerk(this.num_segs, 0.0, 0.0); // +1 seg - velocity ??
        this.num_segs = this.add_segment_const_jerk(this.num_segs, 0.0, 0.0); // +1 seg - velocity ??
        this.num_segs = this.add_segment_const_jerk(this.num_segs, 0.0, 0.0); // +1 seg - velocity ??
        this.num_segs = this.add_segment_const_jerk(this.num_segs, 0.0, 0.0); // +1 seg (15) - velocity ??

        // Calc time at the end of the constant velocity phase
        const t15 = Math.max(0.0, (L - 2.0 * this.segment[this.SEG_SPEED_CHANGE_END].end_pos) / this.segment[this.SEG_SPEED_CHANGE_END].end_vel);

        // Constant velocity phase
        this.num_segs = this.add_segment_const_jerk(this.num_segs, t15, 0.0); // +1 seg (16) - constant velocity

        // Decceleration s-curve phase
        this.num_segs = this.add_segments_jerk(this.num_segs, tj, -Jm, t6);  // +3 segs (19) - decreasing acceleration phase (negative jerk)
        this.num_segs = this.add_segment_const_jerk(this.num_segs, t4, 0.0); // +1 seg (20) - constant acceleration (zero jerk)
        this.num_segs = this.add_segments_jerk(this.num_segs, tj, Jm, t2);   // +3 segs (23) - increasing acceleration phase (positive jerk)

        // remove numerical errors
        this.segment[this.SEG_DECEL_END].end_accel = 0.0;
        this.segment[this.SEG_DECEL_END].end_vel = 0.0;
    }

    // generate three consecutive segments forming a jerk profile
    // the index variable is the position within the path array that this jerk profile should be added
    // the index is incremented to reference the next segment in the array after the jerk profile
    add_segments_jerk(index, tj, Jm, Tcj)
    {
        index = this.add_segment_incr_jerk(index, tj, Jm);
        index = this.add_segment_const_jerk(index, Tcj, Jm);
        index = this.add_segment_decr_jerk(index, tj, Jm);
        return index
    }

    // generate constant jerk time segment
    // calculate the information needed to populate the constant jerk segment from the segment duration tj and jerk J0
    // the index variable is the position of this segment in the path array and is incremented to reference the next segment in the array
    add_segment_const_jerk(index, tj, J0)
    {
        // if no time increase copy previous segment, i.e. there is no constant jerk segment we are going straight
        // into an s-curve jerk segment so we can just copy the exit conditions of the last segment which means we will essentially skip over it at run time
        if (!is_positive(tj)) {
            index = this.add_segment(index,
                                    this.segment[index - 1].end_time,
                                    SegmentType.CONSTANT_JERK,
                                    J0,
                                    this.segment[index - 1].end_accel,
                                    this.segment[index - 1].end_vel,
                                    this.segment[index - 1].end_pos);
            return index;
        }

        // Calculate the segment exit conditions
        const J = J0;
        const T = this.segment[index - 1].end_time + tj;
        const A = this.segment[index - 1].end_accel + J0 * tj;
        const V = this.segment[index - 1].end_vel + this.segment[index - 1].end_accel * tj + 0.5 * J0 * tj**2;
        const P = this.segment[index - 1].end_pos + this.segment[index - 1].end_vel * tj + 0.5 * this.segment[index - 1].end_accel * tj**2 + (1.0 / 6.0) * J0 * Math.pow(tj, 3.0);
        index = this.add_segment(index, T, SegmentType.CONSTANT_JERK, J, A, V, P);
        return index;
    }

    // generate increasing jerk magnitude time segment based on a raised cosine profile
    // calculate the information needed to populate the increasing jerk magnitude segment from the segment duration tj and jerk magnitude Jm
    // the index variable is the position of this segment in the path array and is incremented to reference the next segment in the array
    add_segment_incr_jerk(index, tj, Jm)
    {
        // if no time increase copy previous segment
        if (!is_positive(tj)) {
            index = this.add_segment(index,
                        this.segment[index - 1].end_time,
                        SegmentType.CONSTANT_JERK,
                        0.0,
                        this.segment[index - 1].end_accel,
                        this.segment[index - 1].end_vel,
                        this.segment[index - 1].end_pos);
            return index;
        }
        const Beta = M_PI / tj;
        const Alpha = Jm * 0.5;
        const AT = Alpha * tj;
        const VT = Alpha * (tj * 0.5 - 2.0 / Beta**2);
        const PT = Alpha * ((-1.0 / Beta**2) * tj + (1.0 / 6.0) * Math.pow(tj, 3.0));

        const J = Jm;
        const T = this.segment[index - 1].end_time + tj;
        const A = this.segment[index - 1].end_accel + AT;
        const V = this.segment[index - 1].end_vel + this.segment[index - 1].end_accel * tj + VT;
        const P = this.segment[index - 1].end_pos + this.segment[index - 1].end_vel * tj + 0.5 * this.segment[index - 1].end_accel * tj**2 + PT;
        index = this.add_segment(index, T, SegmentType.POSITIVE_JERK, J, A, V, P);
        return index;
    }

    // generate decreasing jerk magnitude time segment based on a raised cosine profile
    // calculate the information needed to populate the decreasing jerk magnitude segment from the segment duration tj and jerk magnitude Jm
    // the index variable is the position of this segment in the path and is incremented to reference the next segment in the array
    add_segment_decr_jerk(index, tj, Jm)
    {
        // if no time increase copy previous segment
        if (!is_positive(tj)) {
            index = this.add_segment(index,
                        this.segment[index - 1].end_time,
                        SegmentType.CONSTANT_JERK,
                        0.0,
                        this.segment[index - 1].end_accel,
                        this.segment[index - 1].end_vel,
                        this.segment[index - 1].end_pos);
            return index;
        }

        const Beta = M_PI / tj;
        const Alpha = Jm * 0.5;
        const AT = Alpha * tj;
        const VT = Alpha * (tj**2 * 0.5 - 2.0 / Beta**2);
        const PT = Alpha * ((-1.0 / Beta**2) * tj + (1.0 / 6.0) * Math.pow(tj, 3.0));
        const A2T = Jm * tj;
        const V2T = Jm * tj**2;
        const P2T = Alpha * ((-1.0 / Beta**2) * 2.0 * tj + (4.0 / 3.0) * Math.pow(tj, 3.0));

        const J = Jm;
        const T = this.segment[index - 1].end_time + tj;
        const A = (this.segment[index - 1].end_accel - AT) + A2T;
        const V = (this.segment[index - 1].end_vel - VT) + (this.segment[index - 1].end_accel - AT) * tj + V2T;
        const P = (this.segment[index - 1].end_pos - PT) + (this.segment[index - 1].end_vel - VT) * tj + 0.5 * (this.segment[index - 1].end_accel - AT) * tj**2 + P2T;
        index = this.add_segment(index, T, SegmentType.NEGATIVE_JERK, J, A, V, P);
        return index;
    }

    // calculate the segment times for the trigonometric S-Curve path defined by:
    // Sm - duration of the raised cosine jerk profile
    // Jm - maximum value of the raised cosine jerk profile
    // V0 - initial velocity magnitude
    // Am - maximum constant acceleration
    // Vm - maximum constant velocity
    // L - Length of the path
    // tj_out, t2_out, t4_out, t6_out are the segment durations needed to achieve the kinematic path specified by the input variables

    // TODO: NEED to understand this function better.  Will circle back after I have understood how the segments are constructed.
    calculate_path(Sm, Jm, V0, Am, Vm, L)
    {
        // init outputs
        let Jm_out = 0.0;
        let tj_out = 0.0;
        let t2_out = 0.0;
        let t4_out = 0.0;
        let t6_out = 0.0;

        // check for invalid arguments
        if (!is_positive(Sm) || !is_positive(Jm) || !is_positive(Am) || !is_positive(Vm) || !is_positive(L)) {
            throw new Error("SCurve::calculate_path invalid inputs\n");
        }

        if (V0 >= Vm) {
            // no velocity change so all segments as zero length ??????
            return [Jm_out, tj_out, t2_out, t4_out, t6_out];
        }

        // We have been given the maximum snap value (sm).  Due to the known kinematic profile (raised cosine curve)
        // Sm is known to occure at time (t) = tj / 2.  Putting this time into the snap equation to obtain an expression for the
        // the max snap, and rearranging to find the time periosd tj:
        let tj = Jm * M_PI / (2 * Sm);

        const accel_1 = (Vm - V0) / (2.0 * tj);
        const accel_2 = (L + 4.0 * V0 * tj) / (4.0 * tj**2);
        let At = Math.min( Math.min(Am, accel_1), accel_2);

        if (Math.abs(At) < Jm * tj) {
            if (is_zero(V0)) {
                // we do not have a solution for non-zero initial velocity
                let min_tj = Math.min( tj, Math.pow((L * M_PI) / (8.0 * Sm), 1.0/4.0) );
                min_tj = Math.min( min_tj , Math.pow((Vm * M_PI) / (4.0 * Sm), 1.0/3.0) );
                min_tj = Math.min( min_tj, Math.sqrt((Am * M_PI) / (2.0 * Sm)) );
                tj = min_tj;
                Jm = 2.0 * Sm * tj / M_PI;
                Am = Jm * tj;

            } else {
                // When doing speed change we use fixed tj and adjust Jm for small changes
                Am = At;
                Jm = Am / tj;
            }

            const exceed_amax = (Vm <= V0 + 2.0 * Am * tj); // ?????
            const exceed_vmax = (L <= 4.0 * V0 * tj + 4.0 * Am * tj**2); // ??????
            if (exceed_vmax || exceed_vmax) {
                // solution = 0 - t6 t4 t2 = 0 0 0
                t2_out = 0.0;
                t4_out = 0.0;
                t6_out = 0.0;

            } else {
                // solution = 2 - t6 t4 t2 = 0 1 0
                t2_out = 0.0;
                t4_out = Math.min(-(V0 - Vm + Am * tj + (Am * Am) / Jm) / Am, Math.max(((Am * Am) * (-3.0 / 2.0) + Math.sqrt((Am * Am * Am * Am) * (1.0 / 4.0) + (Jm * Jm) * (V0 * V0) + (Am * Am) * (Jm * Jm) * (tj * tj) * (1.0 / 4.0) + Am * (Jm * Jm) * L * 2.0 - (Am * Am) * Jm * V0 + (Am * Am * Am) * Jm * tj * (1.0 / 2.0) - Am * (Jm * Jm) * V0 * tj) - Jm * V0 - Am * Jm * tj * (3.0 / 2.0)) / (Am * Jm), ((Am * Am) * (-3.0 / 2.0) - Math.sqrt((Am * Am * Am * Am) * (1.0 / 4.0) + (Jm * Jm) * (V0 * V0) + (Am * Am) * (Jm * Jm) * (tj * tj) * (1.0 / 4.0) + Am * (Jm * Jm) * L * 2.0 - (Am * Am) * Jm * V0 + (Am * Am * Am) * Jm * tj * (1.0 / 2.0) - Am * (Jm * Jm) * V0 * tj) - Jm * V0 - Am * Jm * tj * (3.0 / 2.0)) / (Am * Jm)));
                t4_out = Math.max(t4_out, 0.0);
                t6_out = 0.0;
            }

        } else {
            if ((Vm < V0 + Am * tj + (Am * Am) / Jm) || (L < 1.0 / (Jm * Jm) * (Am * Am * Am + Am * Jm * (V0 * 2.0 + Am * tj * 2.0)) + V0 * tj * 2.0 + Am * (tj * tj))) {
                // solution = 5 - t6 t4 t2 = 1 0 1
                Am = Math.min(Math.min(Am, Math.max(Jm * (tj + Math.sqrt((V0 * -4.0 + Vm * 4.0 + Jm * (tj * tj)) / Jm)) * (-1.0 / 2.0), Jm * (tj - Math.sqrt((V0 * -4.0 + Vm * 4.0 + Jm * (tj * tj)) / Jm)) * (-1.0 / 2.0))), Jm * tj * (-2.0 / 3.0) + ((Jm * Jm) * (tj * tj) * (1.0 / 9.0) - Jm * V0 * (2.0 / 3.0)) * 1.0 / Math.pow(Math.sqrt(Math.pow(- (Jm * Jm) * L * (1.0 / 2.0) + (Jm * Jm * Jm) * (tj * tj * tj) * (8.0 / 2.7E1) - Jm * tj * ((Jm * Jm) * (tj * tj) + Jm * V0 * 2.0) * (1.0 / 3.0) + (Jm * Jm) * V0 * tj, 2.0) - Math.pow((Jm * Jm) * (tj * tj) * (1.0 / 9.0) - Jm * V0 * (2.0 / 3.0), 3.0)) + (Jm * Jm) * L * (1.0 / 2.0) - (Jm * Jm * Jm) * (tj * tj * tj) * (8.0 / 2.7E1) + Jm * tj * ((Jm * Jm) * (tj * tj) + Jm * V0 * 2.0) * (1.0 / 3.0) - (Jm * Jm) * V0 * tj, 1.0 / 3.0) + Math.pow(Math.sqrt(Math.pow(- (Jm * Jm) * L * (1.0 / 2.0) + (Jm * Jm * Jm) * (tj * tj * tj) * (8.0 / 2.7E1) - Jm * tj * ((Jm * Jm) * (tj * tj) + Jm * V0 * 2.0) * (1.0 / 3.0) + (Jm * Jm) * V0 * tj, 2.0) - Math.pow((Jm * Jm) * (tj * tj) * (1.0 / 9.0) - Jm * V0 * (2.0 / 3.0), 3.0)) + (Jm * Jm) * L * (1.0 / 2.0) - (Jm * Jm * Jm) * (tj * tj * tj) * (8.0 / 2.7E1) + Jm * tj * ((Jm * Jm) * (tj * tj) + Jm * V0 * 2.0) * (1.0 / 3.0) - (Jm * Jm) * V0 * tj, 1.0 / 3.0));
                t2_out = Am / Jm - tj;
                t4_out = 0.0;
                t6_out = t2_out;
            } else {
                // solution = 7 - t6 t4 t2 = 1 1 1
                t2_out = Am / Jm - tj;
                t4_out = Math.min(-(V0 - Vm + Am * tj + (Am * Am) / Jm) / Am, Math.max(((Am * Am) * (-3.0 / 2.0) + Math.sqrt((Am * Am * Am * Am) * (1.0 / 4.0) + (Jm * Jm) * (V0 * V0) + (Am * Am) * (Jm * Jm) * (tj * tj) * (1.0 / 4.0) + Am * (Jm * Jm) * L * 2.0 - (Am * Am) * Jm * V0 + (Am * Am * Am) * Jm * tj * (1.0 / 2.0) - Am * (Jm * Jm) * V0 * tj) - Jm * V0 - Am * Jm * tj * (3.0 / 2.0)) / (Am * Jm), ((Am * Am) * (-3.0 / 2.0) - Math.sqrt((Am * Am * Am * Am) * (1.0 / 4.0) + (Jm * Jm) * (V0 * V0) + (Am * Am) * (Jm * Jm) * (tj * tj) * (1.0 / 4.0) + Am * (Jm * Jm) * L * 2.0 - (Am * Am) * Jm * V0 + (Am * Am * Am) * Jm * tj * (1.0 / 2.0) - Am * (Jm * Jm) * V0 * tj) - Jm * V0 - Am * Jm * tj * (3.0 / 2.0)) / (Am * Jm)));
                t4_out = Math.max(t4_out, 0.0);
                t6_out = t2_out;
            }
        }
        tj_out = tj;
        Jm_out = Jm;

        // check outputs and reset back to zero if necessary
        if (!Number.isFinite(Jm_out) || is_negative(Jm_out) ||
            !Number.isFinite(tj_out) || is_negative(tj_out) ||
            !Number.isFinite(t2_out) || is_negative(t2_out) ||
            !Number.isFinite(t4_out) || is_negative(t4_out) ||
            !Number.isFinite(t6_out) || is_negative(t6_out)) {

            throw new Error("SCurve::calculate_path invalid outputs\n");
        }

        return [Jm_out, tj_out, t2_out, t4_out, t6_out];
    }

    // set the maximum vehicle speed at the origin
    // returns the expected speed at the origin which will always be equal or lower than speed
    set_origin_speed_max(speed)
    {
        // if path is zero length then start speed must be zero
        if (this.num_segs != this.segments_max) {
            return 0.0;
        }

        // avoid re-calculating if unnecessary
        if (is_equal(this.segment[this.SEG_INIT].end_vel, speed)) {
            return speed;
        }

        const Vm = this.segment[this.SEG_ACCEL_END].end_vel;
        const track_length = this.track.length();
        speed = MIN(speed, Vm);

        let [Jm, tj, t2, t4, t6] = this.calculate_path(this.snap_max, this.jerk_max, speed, this.accel_max, Vm, track_length * 0.5);

        // Overwrite the first 8 segements to go from origin speed to constant velocity phase of the s-curve segments
        let seg = this.SEG_INIT;
        seg = this.add_segment(seg, 0.0, SegmentType.CONSTANT_JERK, 0.0, 0.0, speed, 0.0);
        seg = this.add_segments_jerk(seg, tj, Jm, t2);
        seg = this.add_segment_const_jerk(seg, t4, 0.0);
        seg = this.add_segments_jerk(seg, tj, -Jm, t6);

        // remove numerical errors
        this.segment[this.SEG_ACCEL_END].end_accel = 0.0;

        // offset acceleration segment if we can't fit it all into half the original length
        const dPstart = Math.min(0.0, track_length * 0.5 - this.segment[this.SEG_ACCEL_END].end_pos);  // MANNA TODO: I think this maybe backwards
        const dt =  dPstart / this.segment[this.SEG_ACCEL_END].end_vel;
        for (let i = this.SEG_INIT; i <= this.SEG_ACCEL_END; i++) {
            segment[i].end_time += dt;
            segment[i].end_pos += dPstart;
        }

        // add empty speed change segments and constant speed segment  -- TODO: Why do this manually? we have a add_segment function
        for (let i = this.SEG_ACCEL_END+1; i <= this.SEG_SPEED_CHANGE_END; i++) {
            this.segment[i].seg_type = SegmentType.CONSTANT_JERK;
            this.segment[i].jerk_ref = 0.0;
            this.segment[i].end_time = segment[SEG_ACCEL_END].end_time;
            this.segment[i].end_accel = 0.0;
            this.segment[i].end_vel = segment[SEG_ACCEL_END].end_vel;
            this.segment[i].end_pos = segment[SEG_ACCEL_END].end_pos;
        }

        seg = this.SEG_CONST; // TODO count up the segments in the for loop above
        seg = this.add_segment_const_jerk(seg, 0.0, 0.0);

        // now recalc the the deceleration path required to end the s-curve at zero vel and accel????
        // TODO: this also doesn't seem right as we are not at zero velocity (V0 = 0.0 below).  At least I don't think we are.....
        [Jm, tj, t2, t4, t6] = this.calculate_path(this.snap_max, this.jerk_max, 0.0, this.accel_max, this.segment[this.SEG_CONST].end_vel, track_length * 0.5);

        seg = this.add_segments_jerk(seg, tj, -Jm, t6);
        seg = this.add_segment_const_jerk(seg, t4, 0.0);
        seg = this.add_segments_jerk(seg, tj, Jm, t2);

        // remove numerical errors
        this.segment[this.SEG_DECEL_END].end_accel = 0.0;
        this.segment[this.SEG_DECEL_END].end_vel = Math.max(0.0, this.segment[this.SEG_DECEL_END].end_vel);

        // add to constant velocity segment to end at the correct position
        const dP = MAX(0.0, track_length - this.segment[this.SEG_DECEL_END].end_pos);
        const t15 =  dP / this.segment[this.SEG_CONST].end_vel;
        for (let i = this.SEG_CONST; i <= this.SEG_DECEL_END; i++) {
            this.segment[i].end_time += t15;
            this.segment[i].end_pos += dP;
        }

        // catch calculation errors
        if (!valid()) {
            throw new Error("SCurve::set_origin_speed_max invalid path\n");
        }

        return speed;
    }

    // move target location along path from origin to destination
    // prev_leg and next_leg are the paths before and after this path
    // wp_radius is max distance from the waypoint at the apex of the turn
    // fast_waypoint should be true if vehicle will not stop at end of this leg
    // dt is the time increment the vehicle will move along the path
    // target_pos should be set to this segment's origin and it will be updated to the current position target
    // target_vel and target_accel are updated with new targets
    // returns true if vehicle has passed the apex of the corner
    advance_target_along_track(prev_leg, next_leg, wp_radius, accel_corner, fast_waypoint, dt, target_pos, target_vel, target_accel)
    {
        [target_pos, target_vel, target_accel] = prev_leg.move_to_pos_vel_accel(dt, target_pos, target_vel, target_accel);
        [target_pos, target_vel, target_accel] = this.move_from_pos_vel_accel(dt, target_pos, target_vel, target_accel);
        let s_finished = this.finished();

        // // check for change of leg on fast waypoint
        // const time_to_destination = this.get_time_remaining();
        // if (fast_waypoint 
        //     && is_zero(next_leg.get_time_elapsed()) 
        //     && (this.get_time_elapsed() >= this.time_turn_out() - next_leg.time_turn_in()) 
        //     && (this.position_sq >= 0.25 * this.track.length_squared())) {

        //     let turn_pos = new Vector();
        //     turn_pos = turn_pos.subtract(this.track);

        //     let turn_vel = new Vector();
        //     let turn_accel = new Vector();
        //     [turn_pos, turn_vel, turn_accel] = this.move_from_time_pos_vel_accel(get_time_elapsed() + time_to_destination * 0.5, turn_pos, turn_vel, turn_accel);

        //     next_leg.move_from_time_pos_vel_accel(time_to_destination * 0.5, turn_pos, turn_vel, turn_accel);
        //     const speed_min = Math.min(this.get_speed_along_track(), next_leg.get_speed_along_track());
        //     if ((this.get_time_remaining() < next_leg.time_end() * 0.5) && (turn_pos.length() < wp_radius) &&
        //         (new Vector(turn_vel.x, turn_vel.y, 0.0).length() < speed_min) &&
        //         (new Vector(turn_accel.x, turn_accel.y, 0.0).length() < accel_corner))
        //         {
        //         next_leg.move_from_pos_vel_accel(dt, target_pos, target_vel, target_accel);
        //     }

        // } else if (!is_zero(next_leg.get_time_elapsed())) {
        //     [target_pos, target_vel, target_accel] = next_leg.move_from_pos_vel_accel(dt, target_pos, target_vel, target_accel);
        //     if (next_leg.get_time_elapsed() >= this.get_time_remaining()) {
        //         s_finished = true;
        //     }
        // }

        return [s_finished, prev_leg, next_leg, target_pos, target_vel, target_accel];
    }

    // increment time pointer and return the position, velocity and acceleration vectors relative to the destination
    move_to_pos_vel_accel(dt, pos, vel, accel)
    {
        if (!(pos instanceof Vector) || !(vel instanceof Vector) || !(accel instanceof Vector)) {
            throw new Error("pos/vel/accel must be Vectors")
        }

        this.advance_time(dt);

        let [scurve_J1, scurve_A1, scurve_V1, scurve_P1] = this.get_jerk_accel_vel_pos_at_time(this.time);

        pos = pos.add(this.delta_unit.scaler_multiply(scurve_P1));
        vel = vel.add(this.delta_unit.scaler_multiply(scurve_V1));
        accel = accel.add(this.delta_unit.scaler_multiply(scurve_A1));
        this.position_sq = scurve_P1**2;

        // update logging
        this.logger.time.push(this.time);
        this.logger.pos.push(scurve_P1*0.01);
        this.logger.vel.push(scurve_V1*0.01);
        this.logger.accel.push(scurve_A1*0.01);
        this.logger.jerk.push(scurve_J1*0.01);

        // change from relative to destination
        pos = pos.subtract(this.track);

        return [pos, vel, accel]
    }

    // increment time pointer and return the position, velocity and acceleration vectors relative to the origin
    move_from_pos_vel_accel(dt, pos, vel, accel)
    {
        if (!(pos instanceof Vector) || !(vel instanceof Vector) || !(accel instanceof Vector)) {
            throw new Error("pos/vel/accel must be Vectors")
        }

        this.advance_time(dt);

        let [scurve_J1, scurve_A1, scurve_V1, scurve_P1] = this.get_jerk_accel_vel_pos_at_time(this.time);

        pos = pos.add(this.delta_unit.scaler_multiply(scurve_P1));
        vel = vel.add(this.delta_unit.scaler_multiply(scurve_V1));
        accel = accel.add(this.delta_unit.scaler_multiply(scurve_A1));
        this.position_sq = scurve_P1**2;

        // update logging
        this.logger.time.push(this.time);
        this.logger.pos.push(scurve_P1*0.01);
        this.logger.vel.push(scurve_V1*0.01);
        this.logger.accel.push(scurve_A1*0.01);
        this.logger.jerk.push(scurve_J1*0.01);

        return [pos, vel, accel];
    }

    // return the position, velocity and acceleration vectors relative to the origin at a specified time along the path
    move_from_time_pos_vel_accel(time_now, pos, vel, accel)
    {
        if (!(pos instanceof Vector) || !(vel instanceof Vector) || !(accel instanceof Vector)) {
            throw new Error("pos/vel/accel must be Vectors")
        }

        let [scurve_J1, scurve_A1, scurve_V1, scurve_P1] = this.get_jerk_accel_vel_pos_at_time(time_now);

        pos = pos.add(this.delta_unit.scaler_multiply(scurve_P1));
        vel = vel.add(this.delta_unit.scaler_multiply(scurve_V1));
        accel = accel.add(this.delta_unit.scaler_multiply(scurve_A1));

        // update logging
        this.logger.time.push(time_now);
        this.logger.pos.push(scurve_P1*0.01);
        this.logger.vel.push(scurve_V1*0.01);
        this.logger.accel.push(scurve_A1*0.01);
        this.logger.jerk.push(scurve_J1*0.01);

        return [pos, vel, accel];
    }

    // calculate the jerk, acceleration, velocity and position at the provided time
    get_jerk_accel_vel_pos_at_time(time_now)
    {
        // start with zeros as function is void and we want to guarantee all outputs are initialised
        let Jt_out = 0.0;
        let At_out = 0.0;
        let Vt_out = 0.0;
        let Pt_out = 0.0;
        if (this.num_segs != this.segments_max) {
            return [Jt_out, At_out, Vt_out, Pt_out];
        }

        let Jtype;
        let pnt = this.num_segs;
        let Jm, tj, T0, A0, V0, P0;

        // find active segment at time_now
        // TODO: we could look to keep an "active segment" instead of having to find the one we want every time
        for (let i = 0; i < this.num_segs; i++) {
            if (time_now < this.segment[this.num_segs - 1 - i].end_time) {
                pnt = this.num_segs - 1 - i;
            }
        }

        if (pnt == 0) {
            // TODO: not sure we need this, we should have filled in the segments correctly in the first place
            Jtype = SegmentType.CONSTANT_JERK;
            Jm = 0.0;
            tj = 0.0;
            T0 = this.segment[pnt].end_time;
            A0 = this.segment[pnt].end_accel;
            V0 = this.segment[pnt].end_vel;
            P0 = this.segment[pnt].end_pos;

        } else if (pnt == this.num_segs) {
            // TODO: not sure we need this, we should have filled in the segments correctly in the first place
            Jtype = SegmentType.CONSTANT_JERK;
            Jm = 0.0;
            tj = 0.0;
            T0 = this.segment[pnt - 1].end_time;
            A0 = this.segment[pnt - 1].end_accel;
            V0 = this.segment[pnt - 1].end_vel;
            P0 = this.segment[pnt - 1].end_pos;

        } else {
            Jtype = this.segment[pnt].seg_type;
            Jm = this.segment[pnt].jerk_ref;
            tj = this.segment[pnt].end_time - this.segment[pnt - 1].end_time;
            T0 = this.segment[pnt - 1].end_time;
            A0 = this.segment[pnt - 1].end_accel;
            V0 = this.segment[pnt - 1].end_vel;
            P0 = this.segment[pnt - 1].end_pos;
        }

        switch (Jtype) {
            case SegmentType.CONSTANT_JERK:
                [Jt_out, At_out, Vt_out, Pt_out] = this.calc_javp_for_segment_const_jerk(time_now - T0, Jm, A0, V0, P0);
                break;
            case SegmentType.POSITIVE_JERK:
                [Jt_out, At_out, Vt_out, Pt_out] = this.calc_javp_for_segment_incr_jerk(time_now - T0, tj, Jm, A0, V0, P0);
                break;
            case SegmentType.NEGATIVE_JERK:
                [Jt_out, At_out, Vt_out, Pt_out] = this.calc_javp_for_segment_decr_jerk(time_now - T0, tj, Jm, A0, V0, P0);
                break;
        }

        // position along the s-curve leg cannot go backwards
        // TODO: Not really sure we need this?
        Pt_out = Math.max(0.0, Pt_out);

        return [Jt_out, At_out, Vt_out, Pt_out];
    }

    // calculate the jerk, acceleration, velocity and position at time time_now when running the constant jerk time segment
    calc_javp_for_segment_const_jerk(time_now, J0, A0, V0, P0)
    {
        let Jt = J0;
        let At = A0 + J0 * time_now;
        let Vt = V0 + A0 * time_now + 0.5 * J0 * (time_now * time_now);
        let Pt = P0 + V0 * time_now + 0.5 * A0 * (time_now * time_now) + (1.0 / 6.0) * J0 * (time_now * time_now * time_now);
        return [Jt, At, Vt, Pt];
    }

    // Calculate the jerk, acceleration, velocity and position at time time_now when running the increasing jerk magnitude time segment based on a raised cosine profile
    calc_javp_for_segment_incr_jerk(time_now, tj, Jm, A0, V0, P0)
    {
        let Jt = 0.0;
        let At = A0;
        let Vt = V0;
        let Pt = P0;

        if (!is_positive(tj)) {
            return [Jt, At, Vt, Pt];
        }

        const Alpha = Jm * 0.5;
        const Beta = M_PI / tj;
        Jt = Alpha * (1.0 - Math.cos(Beta * time_now));
        At = A0 + Alpha * time_now - (Alpha / Beta) * Math.sin(Beta * time_now);
        Vt = V0 + A0 * time_now + (Alpha * 0.5) * (time_now * time_now) + (Alpha / (Beta * Beta)) * Math.cos(Beta * time_now) - Alpha / (Beta * Beta);
        Pt = P0 + V0 * time_now + 0.5 * A0 * (time_now * time_now) + (-Alpha / (Beta * Beta)) * time_now + Alpha * (time_now * time_now * time_now) / 6.0 + (Alpha / (Beta * Beta * Beta)) * Math.sin(Beta * time_now);
        return [Jt, At, Vt, Pt];
    }

    // Calculate the jerk, acceleration, velocity and position at time time_now when running the decreasing jerk magnitude time segment based on a raised cosine profile
    calc_javp_for_segment_decr_jerk(time_now, tj, Jm, A0, V0, P0)
    {
        let Jt = 0.0;
        let At = A0;
        let Vt = V0;
        let Pt = P0;

        if (!is_positive(tj)) {
            return [Jt, At, Vt, Pt];
        }

        const Alpha = Jm * 0.5;
        const Beta = M_PI / tj;
        const AT = Alpha * tj;
        const VT = Alpha * ((tj * tj) * 0.5 - 2.0 / (Beta * Beta));
        const PT = Alpha * ((-1.0 / (Beta * Beta)) * tj + (1.0 / 6.0) * (tj * tj * tj));
        Jt = Alpha * (1.0 - Math.cos(Beta * (time_now + tj)));
        At = (A0 - AT) + Alpha * (time_now + tj) - (Alpha / Beta) * Math.sin(Beta * (time_now + tj));
        Vt = (V0 - VT) + (A0 - AT) * time_now + 0.5 * Alpha * (time_now + tj) * (time_now + tj) + (Alpha / (Beta * Beta)) * Math.cos(Beta * (time_now + tj)) - Alpha / (Beta * Beta);
        Pt = (P0 - PT) + (V0 - VT) * time_now + 0.5 * (A0 - AT) * (time_now * time_now) + (-Alpha / (Beta * Beta)) * (time_now + tj) + (Alpha / 6.0) * (time_now + tj) * (time_now + tj) * (time_now + tj) + (Alpha / (Beta * Beta * Beta)) * Math.sin(Beta * (time_now + tj));
        return [Jt, At, Vt, Pt];
    }

    is_valid()
    {
        // check number of segments
        if (this.num_segs != this.segments_max) {
            return false;
        }

        for (let i = 0; i < this.num_segs; i++) {
            // jerk_ref should be finite (i.e. not NaN or infinity)
            // time, accel, vel and pos should finite and not negative
            if (!Number.isFinite(this.segment[i].jerk_ref) ||
                !Number.isFinite(this.segment[i].end_time) ||
                !Number.isFinite(this.segment[i].end_accel) ||
                !Number.isFinite(this.segment[i].end_vel) || is_negative(this.segment[i].end_vel) ||
                !Number.isFinite(this.segment[i].end_pos)) {
                return false;
            }

            // time and pos should be increasing
            if (i >= 1) {
                if (is_negative(this.segment[i].end_time - this.segment[i-1].end_time) ||
                    is_negative(this.segment[i].end_pos - this.segment[i-1].end_pos)) {
                    return false;
                }
            }
        }

        // last segment should have zero acceleration
        if (!is_zero(this.segment[this.num_segs-1].end_accel)) {
            return false;
        }

        // if we get this far then the curve must be valid
        return true;
    }

    finished()
    {
        return ((this.time >= this.time_end()) || (this.position_sq >= this.track.length_squared()));
    }

    // time at the end of the sequence
    time_end()
    {
        // GUESS: Number of segments is incomplete so we just send zero time
        if (this.num_segs != this.segments_max) {
            return 0.0;
        }
        const SEG_DECEL_END = this.segments_max-1;
        return this.segment[SEG_DECEL_END].end_time;
    }

    // time left before sequence will complete
    get_time_remaining()
    {
        if (this.num_segs != this.segments_max) {
            return 0.0;
        }
        return this.segment[this.SEG_DECEL_END].end_time - this.time;
    }

    // return the current time elapsed
    get_time_elapsed() { return this.time; }

    // get desired maximum speed along track
    get_speed_along_track() { return this.vel_max; }

    // return time offset used to initiate the turn from leg
    time_turn_out()
    {
        if (this.num_segs != this.segments_max) {
            return 0.0;
        }
        return this.segment[this.SEG_TURN_OUT].end_time;
    }

    // return time offset used to initiate the turn onto leg
    time_turn_in()
    {
        if (this.num_segs != this.segments_max) {
            return 0.0;
        }
        return this.segment[this.SEG_TURN_IN].end_time;
    }

}

// Example Usage
// let curve = new SCurve();
// curve.calculateTrack({ x: 0, y: 0, z: 0 }, { x: 10, y: 10, z: 5 }, 5, 2, 2, 1, 1, 0.1, 0.05);
// curve.debug();

class WPNav {

    constructor() {
        this.wp_jerk = parseFloat(document.getElementById("wpnav_max_jerk").value);

        this.accel_corner = parseFloat(document.getElementById("wpnav_accel_c_max").value);
        this.wp_accel_cmss = parseFloat(document.getElementById("wpnav_accel_xy_max").value);
        this.wp_accel_z_cmss = parseFloat(document.getElementById("wpnav_accel_z_max").value);

        this.wp_speed_cms = parseFloat(document.getElementById("wpnav_speed_xy_max").value);
        this.wp_speed_down_cms = parseFloat(document.getElementById("wpnav_speed_z_dn").value);
        this.wp_speed_up_cms = parseFloat(document.getElementById("wpnav_speed_z_up").value);
        this.wp_desired_speed_xy_cms = 0;

        this.wp_radius_cm = parseFloat(document.getElementById("wpnav_radius").value);

        this.scurve_snap = 0; // calculated later in calc_scurve_jerk_and_snap
        this.scurve_jerk = 0; // calculated later in calc_scurve_jerk_and_snap

        // Just init the locations with zeros. They get overwritten before they are used
        this.origin = new Vector();
        this.destination = new Vector();

        this.attitude_control = new Attitude_Control();
        this.pos_control = new Position_Control();

        this.scurve_prev_leg = new SCurve();
        this.scurve_this_leg = new SCurve();
        this.scurve_next_leg = new SCurve();

        this.track_scalar_dt = 1.0;

        this.flags = ({reached_destination: false, fast_waypoint: false, wp_yaw_set: false});
    }


    wp_and_spline_init(speed_cms = 0.0, stopping_point = new Vector()) {

        // sanity check parameters
        // check _wp_accel_cmss is reasonable
        this.scurve_accel_corner = angle_to_accel(this.pos_control.get_lean_angle_max_cd() * 0.01) * 100;
        this.wp_accel_cmss = Math.min(this.wp_accel_cmss, this.scurve_accel_corner);
        const WPNAV_ACCELERATION = 250;
        this.wp_accel_cmss = (this.wp_accel_cmss <= 0) ? WPNAV_ACCELERATION : this.wp_accel_cmss;

        // check _wp_radius_cm is reasonable
        const WPNAV_WP_RADIUS_MIN = 5.0;
        this.wp_radius_cm = Math.max(this.wp_radius_cm, WPNAV_WP_RADIUS_MIN);

        // check _wp_speed
        const WPNAV_WP_SPEED_MIN = 20.0;
        this.wp_speed_cms = Math.max(this.wp_speed_cms, WPNAV_WP_SPEED_MIN);

        // initialize the desired wp speed if not already done
        this.wp_desired_speed_xy_cms = is_positive(speed_cms) ? speed_cms : this.wp_speed_cms;
        this.wp_desired_speed_xy_cms = Math.max(this.wp_desired_speed_xy_cms, WPNAV_WP_SPEED_MIN);

        this.pos_control.set_max_speed_accel_xy(this.wp_desired_speed_xy_cms, this.wp_accel_cmss);
        this.pos_control.set_max_speed_accel_z(-this.get_default_speed_down(), this.wp_speed_up_cms, this.wp_accel_z_cmss);

        // calculate scurve jerk and jerk time
        if (!is_positive(this.wp_jerk)) {
            this.wp_jerk = this.wp_accel_cmss;
        }
        this.calc_scurve_jerk_and_snap();

        this.scurve_prev_leg.init();
        this.scurve_this_leg.init();
        this.scurve_next_leg.init();
        this.track_scalar_dt = 1.0;

        this.flags.reached_destination = true;
        this.flags.fast_waypoint = false;

        // initialise origin and destination to stopping point
        this.origin = stopping_point;
        this.destination = stopping_point;
    }

    update_wpnav() {

        // see if we need to set new wp active?

        // advance the wp along track
    }

    // helper function to calculate scurve jerk and jerk_time values
    // updates _scurve_jerk and _scurve_snap
    calc_scurve_jerk_and_snap() {
        // calculate jerk
        this.scurve_jerk = Math.min(this.attitude_control.get_ang_vel_roll_max_rads() * GRAVITY_MSS, this.attitude_control.get_ang_vel_pitch_max_rads() * GRAVITY_MSS);
        if (is_zero(this.scurve_jerk)) {
            this.scurve_jerk = this.wp_jerk;
        } else {
            this.scurve_jerk = Math.min(this.scurve_jerk, this.wp_jerk);
        }

        // calculate maximum snap
        // Snap (the rate of change of jerk) uses the attitude control input time constant because multicopters
        // lean to accelerate. This means the change in angle is equivalent to the change in acceleration
        this.scurve_snap = (this.scurve_jerk * M_PI) / (2.0 * Math.max(this.attitude_control.get_input_tc(), 0.1));
        const snap = Math.min(this.attitude_control.get_accel_roll_max_radss(), this.attitude_control.get_accel_pitch_max_radss()) * GRAVITY_MSS;
        if (is_positive(snap)) {
            this.scurve_snap = Math.min(this.scurve_snap, snap);
        }
        // reduce maximum snap by a factor of two from what the aircraft is capable of
        this.scurve_snap *= 0.5;
    }

    /**
     * set waypoint as destination
     * @param {{x: number, y: number, z: number}} destination - Destination vector {x, y, z}
     */
    set_wp_destination (destination) {

        this.scurve_prev_leg.init();
        let origin_speed = 0.0;

        // use previous destination as origin
        this.origin = this.destination;

        // In AP there is a buch of conversions and cases depending on alt frame and whether spline wp or not.
        // That stuff drops out in this simplified version.

        // update destination
        this.destination = destination;

        if (this.flags.fast_waypoint && !this.scurve_next_leg.finished()) {
            // We have a fast WP so we can move onto the next leg without having finished the current scurve leg
            this.scurve_this_leg = this.scurve_next_leg;
        } else {
            this.scurve_this_leg.calculate_track(this.origin, this.destination,
                                                this.pos_control.get_max_speed_xy_cms(), this.pos_control.get_max_speed_up_cms(), this.pos_control.get_max_speed_down_cms(),
                                                this.wp_accel_cmss, this.wp_accel_z_cmss,
                                                this.scurve_snap * 100.0, this.scurve_jerk * 100.0);

            if (!is_zero(origin_speed)) {
                // rebuild start of scurve if we have a non-zero origin speed
                this.scurve_this_leg.set_origin_speed_max(origin_speed);
            }
        }

        this.scurve_next_leg.init();
        this.flags.fast_waypoint = false;   // default waypoint back to slow
        this.flags.reached_destination = false;

        return true; // we can't actually fail this as we don't handle terrain tailes in this js version
    }

    /// advance_wp_target_along_track - move target location along track from origin to destination
    advance_wp_target_along_track(dt)
    {
        // target position, velocity and acceleration from straight line or spline calculators
        let target_pos = new Vector();
        let target_vel = new Vector();
        let target_accel = new Vector();

        let s_finished = false;

        // update target position, velocity and acceleration
        target_pos = this.origin;
        [s_finished, this.scurve_prev_leg, this.scurve_next_leg, target_pos, target_vel, target_accel] = this.scurve_this_leg.advance_target_along_track(this.scurve_prev_leg, this.scurve_next_leg, this.wp_radius_cm, this.scurve_accel_corner, this.flags.fast_waypoint, this.track_scalar_dt * dt, target_pos, target_vel, target_accel);


        // // check if we've reached the waypoint
        // if (!this.flags.reached_destination) {
        //     if (s_finished) {
        //         // "fast" waypoints are complete once the intermediate point reaches the destination
        //         if (this.flags.fast_waypoint) {
        //             this.flags.reached_destination = true;
        //         } else {
        //             // regular waypoints also require the copter to be within the waypoint radius
        //             const dist_to_dest = curr_pos - _destination;
        //             if (dist_to_dest.length_squared() <= sq(_wp_radius_cm)) {
        //                 _flags.reached_destination = true;
        //             }
        //         }
        //     }
        // }

        // successfully advanced along track
        return [target_pos, target_vel, target_accel];
    }

    set_wp_destination_loc (destination) {
        // In AP some conversion is done here to go from location to vector from EKF origin

        return this.set_wp_destination(destination)
    }

    get_default_speed_down() { return Math.abs(this.wp_speed_down_cms); }


}

// An essentially static class that is used to keep the function/class calls looking the same between this JS and the AP C++
class Position_Control {
    constructor()
    {
        this.attitude_control = new Attitude_Control();
        this.pid_accel_z = new PSC_PID();

        this.jerk_max_xy_cmsss = 0;
        this.shaping_jerk_xy = parseFloat(document.getElementById("psc_max_jerk_xy").value);
        this.jerk_max_z_cmsss = 0;
        this.shaping_jerk_z = parseFloat(document.getElementById("psc_max_jerk_z").value);

        this.accel_max_xy_cmss = 0;
        this.accel_max_z_cmss = 0;

        this.vel_max_xy_cms = 0;
        this.vel_max_down_cms = 0;
        this.vel_max_up_cms = 0;

        this.lean_angle_max = parseFloat(document.getElementById("psc_lean_angle_max").value);
    }

    get_lean_angle_max_cd()
    {
        if (!is_positive(this.lean_angle_max)) {
            return this.attitude_control.lean_angle_max_cd();
        }
        return lean_angle_max * 100.0;
    }

    get_max_speed_xy_cms() { return this.vel_max_xy_cms; }
    get_max_speed_up_cms() { return this.vel_max_up_cms; }
    get_max_speed_down_cms() { return this.vel_max_down_cms; }

    /// set_max_speed_accel_xy - set the maximum horizontal speed in cm/s and acceleration in cm/s/s
    ///     This function only needs to be called if using the kinematic shaping.
    ///     This can be done at any time as changes in these parameters are handled smoothly
    ///     by the kinematic shaping.
    set_max_speed_accel_xy(speed_cms, accel_cmss)
    {
        this.vel_max_xy_cms = speed_cms;
        this.accel_max_xy_cmss = accel_cmss;

        // ensure the horizontal jerk is less than the vehicle is capable of
        const jerk_max_cmsss = Math.min(this.attitude_control.get_ang_vel_roll_max_rads(), this.attitude_control.get_ang_vel_pitch_max_rads()) * GRAVITY_MSS * 100.0;
        const snap_max_cmssss = Math.min(this.attitude_control.get_accel_roll_max_radss(), this.attitude_control.get_accel_pitch_max_radss()) * GRAVITY_MSS * 100.0;

        // get specified jerk limit
        this.jerk_max_xy_cmsss = this.shaping_jerk_xy * 100.0;

        // limit maximum jerk based on maximum angular rate
        if (is_positive(jerk_max_cmsss) && this.attitude_control.get_bf_feedforward()) {
            this.jerk_max_xy_cmsss = Math.min(this.jerk_max_xy_cmsss, jerk_max_cmsss);
        }

        // limit maximum jerk to maximum possible average jerk based on angular acceleration
        if (is_positive(snap_max_cmssss) && this.attitude_control.get_bf_feedforward()) {
            this.jerk_max_xy_cmsss = Math.min(0.5 * Math.sqrt(this.accel_max_xy_cmss * snap_max_cmssss), this.jerk_max_xy_cmsss);
        }
    }

    /// set_max_speed_accel_z - set the maximum vertical speed in cm/s and acceleration in cm/s/s
    ///     speed_down can be positive or negative but will always be interpreted as a descent speed.
    ///     This function only needs to be called if using the kinematic shaping.
    ///     This can be done at any time as changes in these parameters are handled smoothly
    ///     by the kinematic shaping.
    set_max_speed_accel_z(speed_down, speed_up, accel_cmss)
    {
        // ensure speed_down is always negative
        speed_down = -Math.abs(speed_down);

        // sanity check and update
        if (is_negative(speed_down)) {
            this.vel_max_down_cms = speed_down;
        }
        if (is_positive(speed_up)) {
            this.vel_max_up_cms = speed_up;
        }
        if (is_positive(accel_cmss)) {
            this.accel_max_z_cmss = accel_cmss;
        }

        // ensure the vertical Jerk is not limited by the filters in the Z accel PID object
        this.jerk_max_z_cmsss = this.shaping_jerk_z * 100.0;
        if (is_positive(this.pid_accel_z.filt_T_hz())) {
            this.jerk_max_z_cmsss = Math.min(this.jerk_max_z_cmsss, Math.min(GRAVITY_MSS * 100.0, this.accel_max_z_cmss) * (M_2PI * this.pid_accel_z.filt_T_hz()) / 5.0);
        }
        if (is_positive(this.pid_accel_z.filt_E_hz())) {
            this.jerk_max_z_cmsss = Math.min(this.jerk_max_z_cmsss, Math.min(GRAVITY_MSS * 100.0, this.accel_max_z_cmss) * (M_2PI * this.pid_accel_z.filt_E_hz()) / 5.0);
        }
    }

}

// An essentially static class that is used to keep the function/class calls looking the same between this JS and the AP C++
class Attitude_Control {
    constructor()
    {
        this.accel_roll_max = parseFloat(document.getElementById("atc_accel_roll_max").value);
        this.accel_pitch_max = parseFloat(document.getElementById("atc_accel_pitch_max").value);
        this.ang_vel_roll_max = parseFloat(document.getElementById("atc_roll_rate_max").value);
        this.ang_vel_pitch_max = parseFloat(document.getElementById("atc_pitch_rate_max").value);
        this.angle_max = parseFloat(document.getElementById("lean_angle_max").value);
        this.input_tc = parseFloat(document.getElementById("atc_input_tc").value);
        this.rate_bf_ff_enabled = parseFloat(document.getElementById("atc_rate_ff_enable").value) > 0;
    }

    lean_angle_max_cd() { return this.angle_max; }

    get_input_tc() { return this.input_tc; }

    get_ang_vel_roll_max_rads() { return radians(this.ang_vel_roll_max); }

    get_ang_vel_pitch_max_rads() { return radians(this.ang_vel_pitch_max); }

    get_accel_roll_max_radss() { return radians(this.accel_roll_max * 0.01); }

    get_accel_pitch_max_radss() { return radians(this.accel_pitch_max * 0.01); }

    get_bf_feedforward() { return this.rate_bf_ff_enabled; }

}

// An essentially static class that is used to keep the function/class calls looking the same between this JS and the AP C++
class PSC_PID {
    constructor() {
        this.filt_targ_hz = parseFloat(document.getElementById("psc_target_filter_cutoff").value);
        this.filt_err_hz = parseFloat(document.getElementById("psc_error_filter_cutoff").value);
    }

    filt_T_hz() { return this.filt_targ_hz; }
    filt_E_hz() { return this.filt_err_hz; }
}



class Vector {
    constructor(x = 0.0, y = 0.0, z = 0.0) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    get_array() {
        return [this.x, this.y, this.z];
    }

    zero() {
        this.x = 0.0;
        this.y = 0.0;
        this.z = 0.0;
    }

    is_zero() {
        return is_zero(this.x) && is_zero(this.y) && is_zero(this.z);
    }

    length_squared() {
        return this.x * this.x + this.y * this.y + this.z * this.z;
    }

    length() {
        return Math.sqrt(this.length_squared());
    }

    add(v2) {
        return new Vector(this.x + v2.x,
                          this.y + v2.y,
                          this.z + v2.z);
    }

    subtract(v2) {
        return new Vector(this.x - v2.x,
                          this.y - v2.y,
                          this.z - v2.z);
    }

    scaler_multiply(s) {
        return new Vector(this.x * s,
                          this.y * s,
                          this.z * s);
    }

    // scaler_divide(s) {
    //     return new Vector(this.x / s,
    //                       this.y / s,
    //                       this.z / s);
    // }

    normalize() {
        const v_length = this.length();
        return new Vector(this.x / v_length,
                          this.y / v_length,
                          this.z / v_length);
    }
}


function update()
{

    let point1 = new Vector(parseFloat(document.getElementById("first_wp_x").value),
                            parseFloat(document.getElementById("first_wp_y").value),
                            parseFloat(document.getElementById("first_wp_z").value));

    let point2 = new Vector(parseFloat(document.getElementById("curr_wp_x").value),
                            parseFloat(document.getElementById("curr_wp_y").value),
                            parseFloat(document.getElementById("curr_wp_z").value));

    let point3 = new Vector(parseFloat(document.getElementById("next_wp_x").value),
                            parseFloat(document.getElementById("next_wp_y").value),
                            parseFloat(document.getElementById("next_wp_z").value));

    let point4 = new Vector(parseFloat(document.getElementById("last_wp_x").value),
                            parseFloat(document.getElementById("last_wp_y").value),
                            parseFloat(document.getElementById("last_wp_z").value));



    // construct "waypoint nav", giving it a starting location.
    let wp_nav = new WPNav()

    // wp_start()
    // This pretty much follows the mode auto logic in copter
    let point1_cm = point1.scaler_multiply(100.0)
    wp_nav.wp_and_spline_init(0.0, point1_cm);

    let point2_cm = point2.scaler_multiply(100.0)
    wp_nav.set_wp_destination_loc(point2_cm);

    // now in wp_run()
    let t = 0.0;
    const dt = 1/100;
    const T = 100;
    const n_steps = Math.floor(T/dt);
    let s_pos = [];
    let s_vel = [];
    let s_accel = [];
    let time = [];
    for (let i = 0; i < n_steps; i++) {
        let [pos_cm, vel_cms, accel_cmss] = wp_nav.advance_wp_target_along_track(dt)


        // logging 3D kinematics to add to the 3D plot
        t += dt;
        s_pos.push(pos_cm.scaler_multiply(0.01));         // (m)
        s_vel.push(vel_cms.scaler_multiply(0.01));        // (m/s)
        s_accel.push(accel_cmss.scaler_multiply(0.01));   // (m/s/s)
        time.push(t);
    }

    // Update plots
    wp_pos_plot.data[0].x = [point1.x, point2.x, point3.x, point4.x];
    wp_pos_plot.data[0].y = [point1.y, point2.y, point3.y, point4.y];
    wp_pos_plot.data[0].z = [point1.z, point2.z, point3.z, point4.z];

    wp_pos_plot.data[1].x = s_pos.map(v => v.x);
    wp_pos_plot.data[1].y = s_pos.map(v => v.y);
    wp_pos_plot.data[1].z = s_pos.map(v => v.z);

    Plotly.redraw("waypoint_plot")

    jerk_plot.data[0].x = wp_nav.scurve_this_leg.logger.time;
    // jerk_plot.data[0].y = wp_nav.scurve_this_leg.logger.jerk.map(v => v.length()); // Total jerk in 1D scurve
    jerk_plot.data[0].y = wp_nav.scurve_this_leg.logger.jerk; // Total jerk in 1D scurve

    Plotly.redraw("jerk_plot")

    accel_plot.data[0].x = wp_nav.scurve_this_leg.logger.time;
    accel_plot.data[0].y = wp_nav.scurve_this_leg.logger.accel; // Total accel in 1D scurve

    Plotly.redraw("accel_plot");

    vel_plot.data[0].x = wp_nav.scurve_this_leg.logger.time;
    vel_plot.data[0].y = wp_nav.scurve_this_leg.logger.vel; // Total vel in 1D scurve

    Plotly.redraw("vel_plot");

    pos_plot.data[0].x = wp_nav.scurve_this_leg.logger.time;
    pos_plot.data[0].y = wp_nav.scurve_this_leg.logger.pos; // Total pos in 1D scurve

    Plotly.redraw("pos_plot");

}


