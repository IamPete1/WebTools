const M_PI = Math.PI
const M_2PI = M_PI * 2.0

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
            xaxis: { title: {text: "East (m)" } },
            yaxis: { title: {text: "North (m)" } },
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
function isZero(value)
{
    return !(Math.abs(value) > 0.0);
}

function isPositive(value)
{
    return value > 0;
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

function isEqual (v1, v2)
{
    return isZero(v1-v2);
}


class SCurve {
    constructor() {
        this.init();
    }

    // Initialization function
    init() {
        this.snap_max = 0.0;
        this.jerk_max = 0.0;
        this.accel_max = 0.0;
        this.vel_max = 0.0;
        this.time = 0.0;
        this.num_segs = 0;
        this.segment = [];
        this.track = { x: 0, y: 0, z: 0 };
        this.delta_unit = { x: 0, y: 0, z: 0 };
        this.position_sq = 0.0;
    }

    // Calculate track motion profile between two 3D points
    calculateTrack(origin, destination, speed_xy, speed_up, speed_down, accel_xy, accel_z, snap_maximum, jerk_maximum) {
        this.init();

        // Compute vector between origin and destination
        let track_temp = subtractVectors(destination, origin);
        if (isZero(track_temp) || isZero(vectorLengthSquared(track_temp))) {
            return;
        }

        this.snap_max = snap_maximum;
        this.jerk_max = jerk_maximum;

        this.setKinematicLimits(origin, destination, speed_xy, speed_up, speed_down, accel_xy, accel_z);

        if (!isPositive(this.snap_max) || !isPositive(this.jerk_max) || 
            !isPositive(this.accel_max) || !isPositive(this.vel_max)) {
            console.error("SCurve: Invalid kinematic parameters");
            return;
        }

        this.track = track_temp;
        const track_length = vectorLength(this.track);

        if (isZero(track_length)) {
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

    // Add trajectory segments
    addSegments(track_length) {
        if (this.isZero(track_length)) return;

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

    // Calculate jerk, acceleration, velocity, position for a segment
    getJAVPAtTime(time_now) {
        let Jt = 0, At = 0, Vt = 0, Pt = 0;

        for (let segment of this.segment) {
            if (time_now < segment.duration) {
                switch (segment.type) {
                    case "CONST_JERK":
                        Jt = segment.jerk;
                        At += Jt * time_now;
                        Vt += At * time_now + 0.5 * Jt * Math.pow(time_now, 2);
                        Pt += Vt * time_now + 0.5 * At * Math.pow(time_now, 2) + (1 / 6) * Jt * Math.pow(time_now, 3);
                        break;

                    case "INCR_JERK":
                        Jt = segment.jerk * (1 - Math.cos(Math.PI * time_now / segment.duration));
                        At += Jt * time_now;
                        Vt += At * time_now;
                        Pt += Vt * time_now;
                        break;

                    case "DECR_JERK":
                        Jt = -segment.jerk * (1 - Math.cos(Math.PI * time_now / segment.duration));
                        At += Jt * time_now;
                        Vt += At * time_now;
                        Pt += Vt * time_now;
                        break;
                }
                break;
            }
            time_now -= segment.duration;
        }

        return { Jt, At, Vt, Pt };
    }

    setOriginSpeedMax(speed) {
        // If the path is zero length, start speed must be zero
        if (this.num_segs !== this.segments_max) {
            return 0.0;
        }
    
        // Avoid recalculating if unnecessary
        if (isEqual(this.segment[this.SEG_INIT].end_vel, speed)) {
            return speed;
        }
    
        const Vm = this.segment[this.SEG_ACCEL_END].end_vel;
        const trackLength = vectorLength(this.track);
        speed = Math.min(speed, Vm);
    
        let Jm, tj, t2, t4, t6;
        this.calculatePath(this.snap_max, this.jerk_max, speed, this.accel_max, Vm, trackLength * 0.5, Jm, tj, t2, t4, t6);
    
        let seg = this.SEG_INIT;
        this.addSegment(seg, 0.0, "CONSTANT_JERK", 0.0, 0.0, speed, 0.0);
        this.addSegmentsJerk(seg, tj, Jm, t2);
        this.addSegmentConstJerk(seg, t4, 0.0);
        this.addSegmentsJerk(seg, tj, -Jm, t6);
    
        // Remove numerical errors
        this.segment[this.SEG_ACCEL_END].end_accel = 0.0;
    
        // Offset acceleration segment if it can't fit in half the original length
        const dPstart = Math.min(0.0, trackLength * 0.5 - this.segment[this.SEG_ACCEL_END].end_pos);
        const dt = dPstart / this.segment[this.SEG_ACCEL_END].end_vel;
        for (let i = this.SEG_INIT; i <= this.SEG_ACCEL_END; i++) {
            this.segment[i].end_time += dt;
            this.segment[i].end_pos += dPstart;
        }
    
        // Add empty speed change segments and constant speed segment
        for (let i = this.SEG_ACCEL_END + 1; i <= this.SEG_SPEED_CHANGE_END; i++) {
            this.segment[i].seg_type = "CONSTANT_JERK";
            this.segment[i].jerk_ref = 0.0;
            this.segment[i].end_time = this.segment[this.SEG_ACCEL_END].end_time;
            this.segment[i].end_accel = 0.0;
            this.segment[i].end_vel = this.segment[this.SEG_ACCEL_END].end_vel;
            this.segment[i].end_pos = this.segment[this.SEG_ACCEL_END].end_pos;
        }
    
        seg = this.SEG_CONST;
        this.addSegmentConstJerk(seg, 0.0, 0.0);
    
        this.calculatePath(this.snap_max, this.jerk_max, 0.0, this.accel_max, this.segment[this.SEG_CONST].end_vel, trackLength * 0.5, Jm, tj, t2, t4, t6);
    
        this.addSegmentsJerk(seg, tj, -Jm, t6);
        this.addSegmentConstJerk(seg, t4, 0.0);
        this.addSegmentsJerk(seg, tj, Jm, t2);
    
        // Remove numerical errors
        this.segment[this.SEG_DECEL_END].end_accel = 0.0;
        this.segment[this.SEG_DECEL_END].end_vel = Math.max(0.0, this.segment[this.SEG_DECEL_END].end_vel);
    
        // Adjust constant velocity segment to end at the correct position
        const dP = Math.max(0.0, trackLength - this.segment[this.SEG_DECEL_END].end_pos);
        const t15 = dP / this.segment[this.SEG_CONST].end_vel;
        for (let i = this.SEG_CONST; i <= this.SEG_DECEL_END; i++) {
            this.segment[i].end_time += t15;
            this.segment[i].end_pos += dP;
        }
    
        // Catch calculation errors
        if (!this.isValid()) {
            console.error("SCurve::setOriginSpeedMax invalid path");
            this.init();
            return 0.0;
        }
    
        return speed;
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

class Vector {
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    get_array() {
        return [this.x, this.y, this.z];
    }
}


function run_flare()
{

    let point1 = new Vector(parseFloat(document.getElementById("last_wp_x").value),
                            parseFloat(document.getElementById("last_wp_y").value),
                            parseFloat(document.getElementById("last_wp_z").value));

    let point2 = new Vector(parseFloat(document.getElementById("curr_wp_x").value),
                            parseFloat(document.getElementById("curr_wp_y").value),
                            parseFloat(document.getElementById("curr_wp_z").value));

    let point3 = new Vector(parseFloat(document.getElementById("next_wp_x").value),
                            parseFloat(document.getElementById("next_wp_y").value),
                            parseFloat(document.getElementById("next_wp_z").value));

    const snap_max = parseFloat(document.getElementById("snap_max").value);
    const jerk_max = parseFloat(document.getElementById("max_jerk").value);
    const accel_z_max = parseFloat(document.getElementById("accel_z_max").value);
    const accel_xy_max = parseFloat(document.getElementById("accel_xy_max").value);
    const speed_xy_max = parseFloat(document.getElementById("speed_xy_max").value);
    const speed_z_up = parseFloat(document.getElementById("speed_z_up").value);
    const speed_z_down = parseFloat(document.getElementById("speed_z_down").value);


    const origin_speed = parseFloat(document.getElementById("initial_vel").value);


    let scurve_this_leg = new SCurve();

    scurve_this_leg.calculateTrack(point1, point2, speed_xy_max, speed_z_up, speed_z_down, accel_xy_max, accel_z_max, snap_max, jerk_max);


    if (!isZero(origin_speed)) {
        // rebuild start of scurve if we have a non-zero origin speed
        scurve_this_leg.setOriginSpeedMax(origin_speed);
    }





    // Update plots
    wp_pos_plot.data[0].x = [point1.x, point2.x, point3.x];
    wp_pos_plot.data[0].y = [point1.y, point2.y, point3.y];
    wp_pos_plot.data[0].z = [point1.z, point2.z, point3.z];

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


