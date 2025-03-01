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
    wp_pos_plot.data = [{ type:'scatter3d',  x:[], y:[], z:[], name: 'path', mode: 'lines+markers', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s³" }]

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

    // Acceleration
    jerk_plot.data = [{ x:[], y:[], name: 'Dumb', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s³" },
                      { x:[], y:[], name: 'Smart', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s³" }]

    jerk_plot.layout = {
        legend: { itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 60, r: 50, t: 20 },
        xaxis: { title: {text: time_scale_label } },
        yaxis: { title: {text: "Jerk (m/s³)" } }
    }

    plot = document.getElementById("jerk_plot")
    Plotly.purge(plot)
    Plotly.newPlot(plot, jerk_plot.data, jerk_plot.layout, { displaylogo: false })

    // Acceleration
    accel_plot.data = [{ x:[], y:[], name: 'Dumb', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s²" },
                       { x:[], y:[], name: 'Smart', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s²" }]

    accel_plot.layout = {
        legend: { itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 60, r: 50, t: 20 },
        xaxis: { title: {text: time_scale_label } },
        yaxis: { title: {text: "Acceleration (m/s²)" } }
    }

    plot = document.getElementById("accel_plot");
    Plotly.purge(plot);
    Plotly.newPlot(plot, accel_plot.data, accel_plot.layout, { displaylogo: false });

    // velocity
    vel_plot.data = [{ x:[], y:[], name: 'Dumb', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s" },
                     { x:[], y:[], name: 'Smart', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m/s" }];

    vel_plot.layout = {
        legend: { itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 60, r: 50, t: 20 },
        xaxis: { title: {text: time_scale_label } },
        yaxis: { title: {text: "Velocity (m/s)" } },
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
    pos_plot.data = [{ x:[], y:[], name: 'Dumb', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m" },
                     { x:[], y:[], name: 'Smart', mode: 'lines', hovertemplate: "<extra></extra>%{x:.2f} s<br>%{y:.2f} m" }]

    pos_plot.layout = {
        legend: { itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 60, r: 50, t: 20 },
        xaxis: { title: {text: time_scale_label } },
        yaxis: { title: {text: "Position (m)" } },
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
    // link_plot_axis_range([
    //     ["jerk_plot", "x", "", jerk_plot],
    //     ["accel_plot", "x", "", accel_plot],
    //     ["vel_plot", "x", "", vel_plot],
    //     ["pos_plot", "x", "", pos_plot],
    // ])

    // Link plot reset
    // link_plot_reset([
    //     ["jerk_plot", jerk_plot],
    //     ["accel_plot", accel_plot],
    //     ["vel_plot", vel_plot],
    //     ["pos_plot", pos_plot],
    // ])
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

function vectorLength(v)
{
    return Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

function vectorLengthSquared(v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

function normalizeVector(v)
{
    let length = this.vectorLength(v);
    return length > 0 ? { x: v.x / length, y: v.y / length, z: v.z / length } : { x: 0, y: 0, z: 0 };
}

function subtractVectors(v1, v2)
{
    return { x: v1.x - v2.x, y: v1.y - v2.y, z: v1.z - v2.z };
}

function isEqual(v1, v2)
{
    return is_zero(v1-v2);
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

function get_lean_angle_max_cd()
{
    const psc_lean_angle_max = parseFloat(document.getElementById("psc_lean_angle_max").value)
    if (psc_lean_angle_max > 0) {
        return psc_lean_angle_max * 100.0;
    }
    return parseFloat(document.getElementById("lean_angle_max").value);
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

class SCurve {
    constructor() {

        this.segments_max = 23; // maximum number of time segments

        // init the array of empty segments
        this.segment = Array.from({ length: this.segments_max }, (_, i) => new Segments());

        this.track = new Vector();
        this.delta_unit = new Vector();

        // this.init();
    }

    // Initialization function
    init() {
        this.snap_max = 0.0;
        this.jerk_max = 0.0;
        this.accel_max = 0.0;
        this.vel_max = 0.0;
        this.time = 0.0;

        this.num_segs = 0;
        this.num_segs = this.add_segment(this.num_segs, 0.0, SegmentType.CONSTANT_JERK, 0.0, 0.0, 0.0, 0.0)

        this.track.zero();
        this.delta_unit.zero();
        this.position_sq = 0.0;
    }

    // increment the internal time
    advance_time(dt) {
        time = MIN(time+dt, time_end());
    }

    // Calculate track motion profile between two 3D points
    calculateTrack(origin, destination, speed_xy, speed_up, speed_down, accel_xy, accel_z, snap_maximum, jerk_maximum) {
        this.init();

        // Compute vector between origin and destination
        let track_temp = subtractVectors(destination, origin);
        if (is_zero(track_temp) || is_zero(vectorLengthSquared(track_temp))) {
            return;
        }

        this.snap_max = snap_maximum;
        this.jerk_max = jerk_maximum;

        this.setKinematicLimits(origin, destination, speed_xy, speed_up, speed_down, accel_xy, accel_z);

        if (!is_positive(this.snap_max) || !is_positive(this.jerk_max) || 
            !is_positive(this.accel_max) || !is_positive(this.vel_max)) {
            console.error("SCurve: Invalid kinematic parameters");
            return;
        }

        this.track = track_temp;
        const track_length = vectorLength(this.track);

        if (is_zero(track_length)) {
            this.delta_unit = { x: 0, y: 0, z: 0 };
        } else {
            this.delta_unit = normalizeVector(this.track);
            this.addSegments(track_length);
        }

        if (!this.isValid()) {
            console.error("SCurve: Invalid path calculation");
            this.init();
        }
    }

    // Set kinematic limits
    setKinematicLimits(origin, destination, speed_xy, speed_up, speed_down, accel_xy, accel_z) {
        speed_xy = Math.abs(speed_xy);
        speed_up = Math.abs(speed_up);
        speed_down = Math.abs(speed_down);
        accel_xy = Math.abs(accel_xy);
        accel_z = Math.abs(accel_z);

        let direction = subtractVectors(destination, origin);
        this.vel_max = this.kinematicLimit(direction, speed_xy, speed_up, speed_down);
        this.accel_max = this.kinematicLimit(direction, accel_xy, accel_z, accel_z);
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

    // Add trajectory segments
    addSegments(track_length) {
        if (this.is_zero(track_length)) return;

        let Jm, tj, t2, t4, t6;
        this.calculatePath(this.snap_max, this.jerk_max, 0.0, this.accel_max, this.vel_max, track_length * 0.5, Jm, tj, t2, t4, t6);

        this.segment.push({ type: "INCR_JERK", jerk: Jm, duration: tj });
        this.segment.push({ type: "CONST_JERK", jerk: 0.0, duration: t4 });
        this.segment.push({ type: "DECR_JERK", jerk: -Jm, duration: t6 });

        // Add a constant speed phase
        this.segment.push({ type: "CONST_VELOCITY", jerk: 0.0, duration: track_length / this.vel_max });

        // Deceleration phase
        this.segment.push({ type: "DECR_JERK", jerk: -Jm, duration: t6 });
        this.segment.push({ type: "CONST_JERK", jerk: 0.0, duration: t4 });
        this.segment.push({ type: "INCR_JERK", jerk: Jm, duration: t2 });
    }



    kinematicLimit(direction, speed_xy, speed_up, speed_down) {
        let xy_speed = Math.sqrt(direction.x * direction.x + direction.y * direction.y);
        let vertical_speed = Math.abs(direction.z) > 0 ? (direction.z > 0 ? speed_up : speed_down) : 0;
        return Math.sqrt(xy_speed * xy_speed + vertical_speed * vertical_speed);
    }

    isValid() {
        return this.segment.length > 0;
    }

    debug() {
        console.log("Segments:", this.segment);
        console.log("Kinematics: ", { snap_max: this.snap_max, jerk_max: this.jerk_max, accel_max: this.accel_max, vel_max: this.vel_max });
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

        this.wp_speed_cms = parseFloat(document.getElementById("wpnav_speed_xy_max").value);

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


    wp_and_spline_init(stopping_point) {

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
        this.leg_dest = destination;

        // set WP destination():
        let scurve_this_leg = new SCurve();

        scurve_this_leg.calculateTrack(leg_origin, point2, this.wp_speed_cms, speed_z_up, speed_z_down, accel_xy_max, accel_z_max, snap_max, jerk_max);

        if (!is_zero(origin_speed)) {
            // rebuild start of scurve if we have a non-zero origin speed
            scurve_this_leg.setOriginSpeedMax(origin_speed);
        }

    }


}

// An essentially static class that is used to keep the function/class calls looking the same between this JS and the AP C++
class Position_Control {
    constructor()
    {
        this.attitude_control = new Attitude_Control();
        this.lean_angle_max = parseFloat(document.getElementById("psc_lean_angle_max").value);
    }

    get_lean_angle_max_cd()
    {
        if (!is_positive(this.lean_angle_max)) {
            return this.attitude_control.lean_angle_max_cd();
        }
        return lean_angle_max * 100.0;
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
    }

    lean_angle_max_cd() { return this.angle_max; }

    get_input_tc() { return this.input_tc; }

    get_ang_vel_roll_max_rads() { return radians(this.ang_vel_roll_max); }

    get_ang_vel_pitch_max_rads() { return radians(this.ang_vel_pitch_max); }

    get_accel_roll_max_radss() { return radians(this.accel_roll_max * 0.01); }

    get_accel_pitch_max_radss() { return radians(this.accel_pitch_max * 0.01); }

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
}


function run_flare()
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

    // const origin_speed = parseFloat(document.getElementById("initial_vel").value);


    // construct "waypoint nav", giving it a starting location.
    let wp_nav = new WPNav()


    // This pretty much follows the mode auto logic in copter
    wp_nav.wp_and_spline_init(point1);






    // Update plots
    wp_pos_plot.data[0].x = [point1.x, point2.x, point3.x, point4.x];
    wp_pos_plot.data[0].y = [point1.y, point2.y, point3.y, point4.y];
    wp_pos_plot.data[0].z = [point1.z, point2.z, point3.z, point4.z];

    Plotly.redraw("waypoint_plot")

    // jerk_plot.data[0].x = t
    // jerk_plot.data[0].y = traj.j
    // jerk_plot.data[1].x = t
    // jerk_plot.data[1].y = calcd_traj.j
    // Plotly.redraw("jerk_plot")

    // accel_plot.data[0].x = t
    // accel_plot.data[0].y = traj.a
    // accel_plot.data[1].x = t
    // accel_plot.data[1].y = calcd_traj.a
    // Plotly.redraw("accel_plot")

    // vel_plot.data[0].x = t
    // vel_plot.data[0].y = traj.v
    // vel_plot.data[1].x = t
    // vel_plot.data[1].y = calcd_traj.v
    // Plotly.redraw("vel_plot")

    // pos_plot.data[0].x = t
    // pos_plot.data[0].y = traj.p
    // pos_plot.data[1].x = t
    // pos_plot.data[1].y = calcd_traj.p
    // Plotly.redraw("pos_plot")

}


