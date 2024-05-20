// MAGFit tool for compass calibration

var DataflashParser
const import_done = import('../JsDataflashParser/parser.js').then((mod) => { DataflashParser = mod.default });

const axis = ['x', 'y', 'z']
const fit_types = { offsets: "Offsets", scale: "Offsets and scale", iron: "Offsets and iron" }
var flight_data = {}
var mag_plot = { x: {}, y:{}, z: {} }
var yaw_change = { mag: {}, att: {} }
var field_length = {}
var motor_comp = {}
var error_plot = {}
var error_bars = {}
const gauss_hovertemplate = "<extra></extra>%{meta}<br>%{x:.2f} s<br>%{y:.2f} mGauss"
function setup_plots() {

    document.title = "ArduPilot MAGFit"

    // Turn off buttons that should not be pressed
    document.getElementById("calculate").disabled = true
    document.getElementById("SaveParams").disabled = true

    const time_scale_label = "Time (s)"

    // Setup flight data plot
    const flight_data_plot = ["Roll", "Pitch", "Throttle", "Altitude"]
    const flight_data_unit = ["deg",  "deg",   "",         "m"]
    flight_data.data = []
    for (let i=0;i<flight_data_plot.length;i++) {
        let axi = "y"
        if (i > 0) {
            axi += (i+1)
        }
        flight_data.data[i] = { mode: "lines",
                                name: flight_data_plot[i],
                                meta: flight_data_plot[i],
                                yaxis: axi,
                                hovertemplate: "<extra></extra>%{meta}<br>%{x:.2f} s<br>%{y:.2f} " + flight_data_unit[i] }
    }

    flight_data.layout = {
        xaxis: { title: {text: time_scale_label },
                 domain: [0.07, 0.93],
                 type: "linear", 
                 zeroline: false, 
                 showline: true, 
                 mirror: true,
                 rangeslider: {} },
        showlegend: false,
        margin: { b: 50, l: 50, r: 50, t: 20 },
    }

    // Set axis to match line colors
    const flight_data_axis_pos = [0, 0.06, 0.94, 1]
    for (let i=0;i<flight_data_plot.length;i++) {
        let axi = "yaxis"
        if (i > 0) {
            axi += (i+1)
        }
        const side = i < 2 ? "left" : "right"
        flight_data.layout[axi] = {title: { text: flight_data_plot[i] },
                                            zeroline: false,
                                            showline: true,
                                            mirror: true,
                                            side: side,
                                            position: flight_data_axis_pos[i],
                                            color: plot_default_color(i) }
        if (i > 0) {
            flight_data.layout[axi].overlaying = 'y'
        }
    }

    var plot = document.getElementById("FlightData")
    Plotly.purge(plot)
    Plotly.newPlot(plot, flight_data.data, flight_data.layout, {displaylogo: false});

    // Update start and end time based on range
    document.getElementById("FlightData").on('plotly_relayout', function(data) {

        function range_update(range) {
            document.getElementById("TimeStart").value = Math.floor(range[0])
            document.getElementById("TimeEnd").value = Math.ceil(range[1])
            if (MAG_Data != null) {
                // If we have data then enable re-calculate on updated range
                set_need_calc(true)
            }
        }

        if ((data['xaxis.range'] !== undefined)) {
            range_update(data['xaxis.range'])
            return
        }

        const range_keys = ['xaxis.range[0]', 'xaxis.range[1]']
        if ((data[range_keys[0]] !== undefined) && (data[range_keys[1]] !== undefined)) {
            range_update([data[range_keys[0]], data[range_keys[1]]])
            return
        }

        const auto_range_key = 'xaxis.autorange'
        if ((data[auto_range_key] !== undefined) && (data[auto_range_key] == true)) {
            range_update([MAG_Data.start_time, MAG_Data.end_time])
        }

    })

    // X, Y, Z component plots
    for (const axi of axis) {
        mag_plot[axi].data = []

        mag_plot[axi].layout = {
            xaxis: {title: {text: time_scale_label }, zeroline: false, showline: true, mirror: true},
            yaxis: {title: {text: "Field " + axi + " (mGauss)" }, zeroline: false, showline: true, mirror: true },
            showlegend: true,
            legend: {itemclick: false, itemdoubleclick: false },
            margin: { b: 50, l: 50, r: 50, t: 20 },
        }

        let plot = document.getElementById("mag_plot_" + axi)
        Plotly.purge(plot)
        Plotly.newPlot(plot, mag_plot[axi].data, mag_plot[axi].layout, {displaylogo: false});
    }

    // Error plot
    // Don't actually ploy any thing in the first index, but it makes the index line up with the above plots
    // That makes the indexes slightly easier, but more importantly makes the auto color line color work the same
    error_plot.data = [ { line: { width: 4, color: "#000000" } }]
    error_plot.layout = {
        xaxis: {title: {text: time_scale_label }, zeroline: false, showline: true, mirror: true},
        yaxis: {title: {text: "Field error (mGauss)" }, zeroline: false, showline: true, mirror: true },
        showlegend: true,
        legend: {itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 50, r: 50, t: 20 },
    }
    plot = document.getElementById("error_plot")
    Plotly.purge(plot)
    Plotly.newPlot(plot, error_plot.data, error_plot.layout, {displaylogo: false});

    // Mag yaw error
    yaw_change.mag.data = [ { line: { width: 4, color: "#000000" } }]
    yaw_change.mag.layout = {
        xaxis: {title: {text: time_scale_label }, zeroline: false, showline: true, mirror: true},
        yaxis: {title: {text: "Change heading<br> New vs existing calibration (deg)" }, zeroline: false, showline: true, mirror: true },
        showlegend: true,
        legend: {itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 50, r: 50, t: 20 },
    }
    plot = document.getElementById("yaw_change_mag")
    Plotly.purge(plot)
    Plotly.newPlot(plot, yaw_change.mag.data, yaw_change.mag.layout, {displaylogo: false});

    // ATT yaw error
    yaw_change.att.data = [ { line: { width: 4, color: "#000000" } }]
    yaw_change.att.layout = {
        xaxis: {title: {text: time_scale_label }, zeroline: false, showline: true, mirror: true},
        yaxis: {title: {text: "Change heading<br> New vs selected attitude source (deg)" }, zeroline: false, showline: true, mirror: true },
        showlegend: true,
        legend: {itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 50, r: 50, t: 20 },
    }
    plot = document.getElementById("yaw_change_att")
    Plotly.purge(plot)
    Plotly.newPlot(plot, yaw_change.att.data, yaw_change.att.layout, {displaylogo: false});

    // Mag field length
    field_length.data = [{
        mode: "lines",
        name: "Expected",
        meta: "Expected",
        line: { width: 4, color: "#000000" },
        hovertemplate: gauss_hovertemplate
    }]
    field_length.layout = {
        xaxis: {title: {text: time_scale_label }, zeroline: false, showline: true, mirror: true},
        yaxis: {title: {text: "Measured field length (mGauss)" }, zeroline: false, showline: true, mirror: true },
        showlegend: true,
        legend: {itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 50, r: 50, t: 20 },
    }
    plot = document.getElementById("field_length")
    Plotly.purge(plot)
    Plotly.newPlot(plot, field_length.data, field_length.layout, {displaylogo: false});

    // Motor compensation source
    motor_comp.data = []
    motor_comp.layout = {
        xaxis: {title: {text: time_scale_label }, zeroline: false, showline: true, mirror: true},
        yaxis: {title: {text: "Current (A)" }, zeroline: false, showline: true, mirror: true },
        showlegend: true,
        legend: {itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 50, r: 50, t: 20 },
    }
    plot = document.getElementById("motor_comp")
    Plotly.purge(plot)
    Plotly.newPlot(plot, motor_comp.data, motor_comp.layout, {displaylogo: false});

    // Error bar graph
    error_bars.data = []
    for (let i=0;i<3;i++) {
        const name = "Mag " + (i + 1)
        error_bars.data[i] = { 
            type: "bar",
            name: name,
            meta: name,
            marker: { color: plot_default_color(i+1) }, 
            hovertemplate: "<extra></extra>%{meta}<br>%{x}<br>%{y:.2f} mGauss"
        }
    }
    error_bars.layout = {
        xaxis: { zeroline: false, showline: true, mirror: true},
        yaxis: { title: {text: "mean field error (mGauss)" }, zeroline: false, showline: true, mirror: true },
        barmode: 'group',
        legend: {itemclick: false, itemdoubleclick: false },
        margin: { b: 50, l: 50, r: 50, t: 20 },
    }

    plot = document.getElementById("error_bars")
    Plotly.purge(plot)
    Plotly.newPlot('error_bars', error_bars.data, error_bars.layout, { modeBarButtonsToRemove: ['lasso2d', 'select2d'], displaylogo: false })

}

const offsets_range = [-1500.0, 1500.0]
const diagonals_range = [0.8, 1.2]
const off_diagonals_range = [-0.2, 0.2]
const scale_range = [0.8, 1.2]
function save_parameters() {

    function param_string(name, value) {
        return name + "," + param_to_string(value) + "\n"
    }
    
    function save_params(names, values) {

        function param_array(names, values) {
            var ret = "";
            for (let i = 0; i < names.length; i++) {
                ret += param_string(names[i], values[i])
            }
            return ret
        }

        var ret = "";
        ret += param_array(names.offsets, values.offsets)
        ret += param_array(names.diagonals, values.diagonals)
        ret += param_array(names.off_diagonals, values.off_diagonals)
        ret += param_array(names.motor, values.motor)

        ret += param_string(names.scale, values.scale)
        ret += param_string(names.orientation, values.orientation)
        return ret
    }

    function check_params(i, names, values, original_values) {

        function check_range(name, value, range) {
            let ret = ""
            if (value > range[1]) {
                ret = name + " " + value + " larger than " + range[1] + "\n"
            } else if (value < range[0]) {
                ret = name + " " + value + " less than " + range[0] + "\n"
            }
            return ret
        }

        function check_array(names, values, range) {
            let ret = ""
            for (let i = 0; i < names.length; i++) {
                ret += check_range(names[i], values[i], range)
            }
            return ret
        }

        let warning = ""

        warning += check_array(names.offsets, values.offsets, offsets_range)
        warning += check_array(names.diagonals, values.diagonals, diagonals_range)
        warning += check_array(names.off_diagonals, values.off_diagonals, off_diagonals_range)
        warning += check_range(names.scale, values.scale, scale_range)

        if (warning != "") {
            warning = "MAG " + (i+1) + " params outside typical range:\n" + warning
        }

        if (original_values.orientation != values.orientation) {
            if (warning != "") {
                warning += "\n"
            }
            warning += "MAG " + (i+1) + " orientation (" + names.orientation + ") changed from " + get_rotation_name(original_values.orientation) + " to " + get_rotation_name(values.orientation) + "\n"
        }

        if (warning == "") {
            return true
        }

        return confirm(warning);
    }

    let params = ""
    let type = 0
    let saved = "Saved:\n"
    for (let i = 0; i < 3; i++) {
        if (MAG_Data[i] == null) {
            continue
        }
        const len = MAG_Data[i].param_selection.length
        for (let j = 0; j < len; j++) {
            if (MAG_Data[i].param_selection[j].show) {
                if (check_params(i, MAG_Data[i].names, MAG_Data[i].param_selection[j], MAG_Data[i].params)) {
                    if (!array_all_equal(MAG_Data[i].param_selection[j].motor, 0.0)) {
                        // Check for conflicting motor compensation types
                        const fit_type = MAG_Data[i].param_selection[j].fit_type
                        if (type == 0) {
                            type = fit_type
                        } else if ((fit_type != 0) && (fit_type != type)) {
                            alert("All compasses must use the same motor fit type, current and throttle compensation cannot be used together")
                            return
                        }
                    }
                    params += save_params(MAG_Data[i].names, MAG_Data[i].param_selection[j])

                    // Set use param if selected
                    const option = document.querySelector("input[name=\"MAG" + i + "use\"]:checked").value
                    if (option != 0) {
                        // Override
                        params += param_string(MAG_Data[i].names.use, option == 1 ? 1 : 0)
                    }

                    saved += "\tCompass " + (i + 1) + ": " + MAG_Data[i].param_selection[j].name + "\n"

                }
                break
            }
        }
    }
    if (params == "") {
        alert("No parameters to save")
        return
    }

    params += param_string("COMPASS_MOTCT", type)

    var blob = new Blob([params], { type: "text/plain;charset=utf-8" });
    saveAs(blob, "MAGFit.param");

    alert(saved)
}

function update_shown_params() {

    for (let i = 0; i < 3; i++) {
        if (MAG_Data[i] == null) {
            continue
        }

        const info_name = "MAG_" + (i + 1) + "_PARAM_INFO"

        function show_params(names, values, fit_name) {
            for (let i = 0; i < 3; i++) {
                parameter_set_value(names.offsets[i], values.offsets[i])
                parameter_set_value(names.diagonals[i], values.diagonals[i])
                parameter_set_value(names.off_diagonals[i], values.off_diagonals[i])
                parameter_set_value(names.motor[i], values.motor[i])
            }
            parameter_set_value(names.scale, values.scale)
            parameter_set_value(names.orientation, values.orientation)
    
            document.getElementById(info_name).replaceChildren(document.createTextNode(fit_name))
        }

        let shown = false

        const len = MAG_Data[i].param_selection.length
        for (let j = 0; j < len; j++) {
            if (MAG_Data[i].param_selection[j].show) {
                show_params(MAG_Data[i].names, MAG_Data[i].param_selection[j], MAG_Data[i].param_selection[j].name)
                shown = true
                break
            }
        }

        if (!shown) {
            // Show existing params
            show_params(MAG_Data[i].names, MAG_Data[i].params, "Existing calibration")
        }
    }

}

function redraw() {


    // Expected field
    for (const axi of axis) {
        mag_plot[axi].data = []
    }
    error_plot.data = [ error_plot.data[0] ]
    yaw_change.mag.data = [ yaw_change.mag.data[0] ]
    yaw_change.att.data = [ yaw_change.att.data[0] ]

    field_length.data = [ field_length.data[0] ]
    field_length.data[0].x = [MAG_Data.start_time, MAG_Data.end_time]
    field_length.data[0].y = [earth_field.intensity * 1000.0, earth_field.intensity * 1000.0]

    for (let i = 0; i < 3; i++) {
        if (MAG_Data[i] == null) {
            continue
        }
        const name = String(i + 1)
        const calibrated_name = name + " calibrated"
        const expected_name = name + " expected"
        MAG_Data[i].param_selection = []

        function setup_plot(data, group_index, fit_name) {
            const show = data.show.checked
            const index = mag_plot.x.data.length / 2
            data.show.setAttribute('data-index', index)

            const group_name = fit_name

            for (const axi of axis) {
                mag_plot[axi].data.push({
                    mode: "lines",
                    name: calibrated_name,
                    meta: calibrated_name,
                    visible: show,
                    legendgroup: group_index,
                    legendgrouptitle: { text: group_name },
                    hovertemplate: gauss_hovertemplate,
                    x: MAG_Data[i].time,
                    y: data[axi]
                })

                mag_plot[axi].data.push({
                    mode: "lines",
                    name: expected_name,
                    meta: expected_name,
                    visible: show,
                    legendgroup: group_index,
                    legendgrouptitle: { text: group_name },
                    hovertemplate: gauss_hovertemplate,
                    x: MAG_Data[i].time,
                    y: data.expected[axi],
                    line: { dash: "dot" },
                })
            }

            error_plot.data.push({
                mode: "lines",
                name: name,
                meta: name,
                visible: show,
                legendgroup: group_index,
                legendgrouptitle: { text: group_name },
                hovertemplate: gauss_hovertemplate,
                x: MAG_Data[i].time,
                y: data.error
            })

            // Calculate yaw change from existing calibration
            // always 0 for existing line, hide plot with 0 line width
            // This means line colors still match
            const existing = (group_index == -1)

            let mag_yaw = 0
            if ("yaw" in data) {
                mag_yaw = array_scale(array_wrap_PI(array_sub(data.yaw, MAG_Data[i].orig.yaw)), 180 / Math.PI)
            }

            const yaw_change_hover = "<extra></extra>%{meta}<br>%{x:.2f} s<br>%{y:.2f} deg"

            yaw_change.mag.data.push({
                mode: "lines",
                name: name,
                meta: name,
                line: { width: existing ? 0 : 2 },
                visible: show,
                showlegend: !existing,
                legendgroup: group_index,
                legendgrouptitle: { text: group_name },
                hoverinfo : existing ? "none" : "all",
                hovertemplate: existing ? "" : yaw_change_hover,
                x: MAG_Data[i].time,
                y: mag_yaw
            })

            let att_yaw = 0
            if ("yaw" in data) {
                att_yaw = array_scale(array_wrap_PI(array_sub(data.yaw, MAG_Data[i].quaternion.yaw)), 180 / Math.PI)
            }

            yaw_change.att.data.push({
                mode: "lines",
                name: name,
                meta: name,
                visible: show,
                legendgroup: group_index,
                legendgrouptitle: { text: group_name },
                hovertemplate: yaw_change_hover,
                x: MAG_Data[i].time,
                y: att_yaw
            })

            // Calculate field length
            let magnitude
            if ("x" in data) {
                const len = data.x.length
                magnitude = new Array(len)
                for (let j = 0; j < len; j++) {
                    magnitude[j] = Math.sqrt(data.x[j]**2 + data.y[j]**2 + data.z[j]**2)
                }
            }

            field_length.data.push({
                mode: "lines",
                name: name,
                meta: name,
                visible: show,
                legendgroup: group_index,
                legendgrouptitle: { text: group_name },
                hovertemplate: gauss_hovertemplate,
                x: MAG_Data[i].time,
                y: magnitude
            })

            if ("params" in data) {
                MAG_Data[i].param_selection.push(Object.assign({
                    // Use same index as plots
                    index: index,
                    name: group_name.replace("<br>", ', '),
                    show
                }, data.params))
            }
        }

        // Existing cal
        setup_plot(MAG_Data[i].orig, -1, "Existing cal")

        for (let j = 0; j < MAG_Data[i].fits.length; j++) {
            const group_index = j*Object.keys(fit_types).length
            setup_plot(MAG_Data[i].fits[j].offsets, group_index + 0, fit_types.offsets + "<br>" + MAG_Data[i].fits[j].name)
            setup_plot(MAG_Data[i].fits[j].scale,   group_index + 1, fit_types.scale + "<br>" +  MAG_Data[i].fits[j].name)
            setup_plot(MAG_Data[i].fits[j].iron,    group_index + 2, fit_types.iron + "<br>" +  MAG_Data[i].fits[j].name)
        }

    }

    const time_range = [ parseFloat(document.getElementById("TimeStart").value),
                         parseFloat(document.getElementById("TimeEnd").value)]

    for (const axi of axis) {
        mag_plot[axi].layout.xaxis.autorange = false
        mag_plot[axi].layout.xaxis.range = time_range
        Plotly.newPlot("mag_plot_" + axi, mag_plot[axi].data, mag_plot[axi].layout, {displaylogo: false});
    }

    error_plot.layout.xaxis.autorange = false
    error_plot.layout.xaxis.range = time_range

    Plotly.newPlot("error_plot", error_plot.data, error_plot.layout, {displaylogo: false});

    yaw_change.mag.layout.xaxis.autorange = false
    yaw_change.mag.layout.xaxis.range = time_range

    Plotly.newPlot("yaw_change_mag", yaw_change.mag.data, yaw_change.mag.layout, {displaylogo: false});

    yaw_change.att.layout.xaxis.autorange = false
    yaw_change.att.layout.xaxis.range = time_range

    Plotly.newPlot("yaw_change_att", yaw_change.att.data, yaw_change.att.layout, {displaylogo: false});

    field_length.layout.xaxis.autorange = false
    field_length.layout.xaxis.range = time_range

    Plotly.newPlot("field_length", field_length.data, field_length.layout, {displaylogo: false});

    motor_comp.layout.xaxis.autorange = false
    motor_comp.layout.xaxis.range = time_range

    Plotly.redraw("motor_comp");

    // Plot error bars
    for (let i = 0; i < 3; i++) {
        if (MAG_Data[i] == null) {
            continue
        }

        error_bars.data[i].y = [MAG_Data[i].orig.mean_error]
        error_bars.data[i].x = ["Existing Calibration"]

        for (let j = 0; j < MAG_Data[i].fits.length; j++) {
            const name =  "<br>" + MAG_Data[i].fits[j].name

            if (MAG_Data[i].fits[j].offsets.valid) {
                error_bars.data[i].x.push("Offsets" + name)
                error_bars.data[i].y.push(MAG_Data[i].fits[j].offsets.mean_error)
            }
            if (MAG_Data[i].fits[j].scale.valid) {
                error_bars.data[i].x.push("Offsets and scale" + name)
                error_bars.data[i].y.push(MAG_Data[i].fits[j].scale.mean_error)
            }
            if (MAG_Data[i].fits[j].iron.valid) {
                error_bars.data[i].x.push("Offsets and iron" + name)
                error_bars.data[i].y.push(MAG_Data[i].fits[j].iron.mean_error)
            }
        }

        error_bars.data[i].visible = true // document.getElementById("MAG" + i + "_SHOW").checked
    }
    Plotly.redraw("error_bars")

    // Clear listeners
    document.getElementById("mag_plot_x").removeAllListeners("plotly_relayout");
    document.getElementById("mag_plot_y").removeAllListeners("plotly_relayout");
    document.getElementById("mag_plot_z").removeAllListeners("plotly_relayout");
    document.getElementById("error_plot").removeAllListeners("plotly_relayout");
    document.getElementById("error_bars").removeAllListeners("plotly_relayout");
    document.getElementById("yaw_change_mag").removeAllListeners("plotly_relayout");
    document.getElementById("yaw_change_att").removeAllListeners("plotly_relayout");
    document.getElementById("field_length").removeAllListeners("plotly_relayout");
    document.getElementById("motor_comp").removeAllListeners("plotly_relayout");

    // Link all time axis
    link_plot_axis_range([
        ["mag_plot_x", "x", "", mag_plot.x],
        ["mag_plot_y", "x", "", mag_plot.y],
        ["mag_plot_z", "x", "", mag_plot.z],
        ["error_plot", "x", "", error_plot],
        ["yaw_change_mag", "x", "", yaw_change.mag],
        ["yaw_change_att", "x", "", yaw_change.att],
        ["field_length", "x", "", field_length],
        ["motor_comp", "x", "", motor_comp],
    ])

    // Link plot reset
    link_plot_reset([
        ["mag_plot_x", mag_plot.x],
        ["mag_plot_y", mag_plot.y],
        ["mag_plot_z", mag_plot.z],
        ["error_plot", error_plot],
        ["error_bars", error_bars],
        ["yaw_change_mag", yaw_change.mag],
        ["yaw_change_att", yaw_change.att],
        ["field_length", field_length],
        ["motor_comp", motor_comp],
    ])

}

function update_hidden(ele) {
    if (ele.dataset.index == null) {
        return
    }
    const index = parseFloat(ele.dataset.index)
    const show = ele.checked

    const axi_index = index * 2
    for (const axi of axis) {
        mag_plot[axi].data[axi_index + 0].visible = show
        mag_plot[axi].data[axi_index + 1].visible = show
        Plotly.redraw("mag_plot_" + axi)
    }

    error_plot.data[index].visible = show
    Plotly.redraw("error_plot")

    yaw_change.mag.data[index].visible = show
    Plotly.redraw("yaw_change_mag")

    yaw_change.att.data[index].visible = show
    Plotly.redraw("yaw_change_att")

    field_length.data[index].visible = show
    Plotly.redraw("field_length")

    // Turn off error bars for any mag that is completely disabled
    for (let i = 0; i < 3; i++) {
        if (MAG_Data[i] == null) {
            continue
        }
        let show_bar = false

        show_bar ||= MAG_Data[i].orig.show.checked
        for (let j = 0; j < MAG_Data[i].fits.length; j++) {
            for (const key of Object.keys(fit_types)) {
                show_bar ||= MAG_Data[i].fits[j][key].show.checked
            }
        }

        error_bars.data[i].visible = show_bar
    }
    Plotly.redraw("error_bars")

    // Add/remove the set param set from list
    for (let i = 0; i < 3; i++) {
        if (MAG_Data[i] == null) {
            continue
        }
        const selection_index = MAG_Data[i].param_selection.findIndex(i => i.index === index)
        if (selection_index == -1) {
            // Not found
            continue
        }

        // set show
        MAG_Data[i].param_selection[selection_index].show = show
        if (show) {
            // Move to first priority
            const selection = MAG_Data[i].param_selection[selection_index]

            // Remove
            MAG_Data[i].param_selection.splice(selection_index, 1)

            // Add to start
            MAG_Data[i].param_selection.splice(0, 0, selection)

        }
    }

    update_shown_params()

}

function scale_valid(scale) {
    const MAX_SCALE_FACTOR = 1.5
    return (scale <= MAX_SCALE_FACTOR) && (scale >= (1/MAX_SCALE_FACTOR))
}

function apply_params(ret, raw, params, motor) {

    // Offsets
    ret.x = array_offset(raw.x, params.offsets[0])
    ret.y = array_offset(raw.y, params.offsets[1])
    ret.z = array_offset(raw.z, params.offsets[2])

    // scale
    if (scale_valid(params.scale)) {
        ret.x = array_scale(ret.x, params.scale)
        ret.y = array_scale(ret.y, params.scale)
        ret.z = array_scale(ret.z, params.scale)
    }

    // Iron
    if (!array_all_equal(params.diagonals, 0.0)) {

        // Vectorized multiplication
        const corrected_x = array_add(array_add( array_scale(ret.x, params.diagonals[0]),     array_scale(ret.y, params.off_diagonals[0])), array_scale(ret.z, params.off_diagonals[1]) )
        const corrected_y = array_add(array_add( array_scale(ret.x, params.off_diagonals[0]), array_scale(ret.y, params.diagonals[1])),     array_scale(ret.z, params.off_diagonals[2]) )
        const corrected_z = array_add(array_add( array_scale(ret.x, params.off_diagonals[1]), array_scale(ret.y, params.off_diagonals[2])), array_scale(ret.z, params.diagonals[2]) )

        ret.x = corrected_x; ret.y = corrected_y; ret.z = corrected_z
    }

    // Motor
    if (motor != null) {
        ret.x = array_add(ret.x, array_scale(motor, params.motor[0]))
        ret.y = array_add(ret.y, array_scale(motor, params.motor[1]))
        ret.z = array_add(ret.z, array_scale(motor, params.motor[2]))
    }
}

// Look through time array and return first index before start time
function find_start_index(time) {
    const start_time = parseFloat(document.getElementById("TimeStart").value)

    var start_index = 0
    for (j = 0; j<time.length; j++) {
        // Move forward start index while time is less than start time
        if (time[j] < start_time) {
            start_index = j
        }
    }
    return start_index
}

// Look through time array and return first index after end time
function find_end_index(time) {
    const end_time = parseFloat(document.getElementById("TimeEnd").value)

    var end_index = 0
    for (j = 0; j<time.length-1; j++) {
        // Move forward end index while time is less than end time
        if (time[j] <= end_time) {
            end_index = j + 1
        }
    }
    return end_index
}

// Run all calculation steps
function calculate() {

    select_body_frame_attitude()

    calculate_bins()

    check_orientation()

    fit()

    redraw()

    update_shown_params()

    set_need_calc(false)

}

// Angle wrap helpers
function wrap_2PI(rad) {
    const PI2 = Math.PI * 2
    return (rad + PI2) % PI2
}

function wrap_PI(rad) {
    let ret = wrap_2PI(rad)
    if (ret > Math.PI) {
        ret -= Math.PI * 2
    }
    return ret
}

function array_wrap_PI(A) {
    const len = A.length
    let ret = new Array(len)
    for (let i = 0; i < len; i++) {
        ret[i] = wrap_PI(A[i])
    }
    return ret
}

function array_wrap_2PI(A) {
    const len = A.length
    let ret = new Array(len)
    for (let i = 0; i < len; i++) {
        ret[i] = wrap_2PI(A[i])
    }
    return ret
}

// Calculate yaw estimate from compass only, tilt correction
function get_yaw(mag_field, tilt_correction) {
    const len = mag_field.x.length
    const declination_rad = earth_field.declination * (Math.PI / 180)

    let yaw = new Array(len)
    for (let i = 0; i < len; i++) {

        // Remote roll and pitch from measured felid (tilt correction)
        const xyz = tilt_correction[i].rotate([mag_field.x[i], mag_field.y[i], mag_field.z[i]])

        // Calculate yaw from measurement
        yaw[i] = wrap_2PI(Math.atan2(-xyz[1],xyz[0]) + declination_rad)

    }

    return yaw
}

// Calculate error weights based on attitude binning
const num_bins = 80
function calculate_bins() {
    const start = performance.now()

    // Fibonacci lattice of unit radius
    const bins = { x: new Array(num_bins), y: new Array(num_bins), z: new Array(num_bins) }
    for (let i = 0; i < num_bins; i++) {
        const k = i + 0.5;

        const phi = Math.acos(1.0 - 2.0 * k / num_bins)
        const theta = Math.PI * (1 + Math.sqrt(5)) * k

        bins.x[i] = Math.cos(theta) * Math.sin(phi)
        bins.y[i] = Math.sin(theta) * Math.sin(phi)
        bins.z[i] = Math.cos(phi)
    }

    for (let i = 0; i < 3; i++) {
        if (MAG_Data[i] == null) {
            continue
        }
        const len = MAG_Data[i].orig.expected.x.length
        MAG_Data[i].bins = new Array(len)

        // Find the closest bin to each point 
        for (let j = 0; j < len; j++) {
            // Convert to unit
            let x = MAG_Data[i].orig.expected.x[j]
            let y = MAG_Data[i].orig.expected.y[j]
            let z = MAG_Data[i].orig.expected.z[j]
            const length = Math.sqrt(x**2 + y**2 + z**2)
            x /= length
            y /= length
            z /= length

            // Check all points
            let min_dist = Infinity
            for (let k = 0; k < num_bins; k++) {
                const dist_sq = (x - bins.x[k])**2 + (y - bins.y[k])**2 + (z - bins.z[k])**2

                if (dist_sq < min_dist) {
                    min_dist = dist_sq
                    MAG_Data[i].bins[j] = k
                }
            }
        }
    }

    const end = performance.now();
    console.log(`Binning took: ${end - start} ms`);
}

function get_weights(bins) {

    const count = new Array(num_bins).fill(0)

    const len = bins.length
    let num_unique_bins = 0
    let total_bins = 0
    for (let i = 0; i < len; i++) {
        if (count[bins[i]] == 0) {
            num_unique_bins++
        }
        count[bins[i]]++
        total_bins++
    }
    const mean_bin_size = total_bins / num_unique_bins
    const coverage = num_unique_bins / num_bins

    let weights = new Array(len).fill(1)

    for (let i = 0; i < len; i++) {
        // Scale by mean_bin_size so that the average weight is 1, this give comparable error magnitude to the un-weighted case
        //weights[i] = mean_bin_size / count[bins[i]]
    }

    return { weights, coverage }
}

function check_orientation() {

    const start = performance.now()

    for (let i = 0; i < 3; i++) {
        if (MAG_Data[i] == null || !MAG_Data[i].rotate) {
            continue
        }

        const option = document.querySelector("input[name=\"MAG" + i + "orientation\"]:checked").value
        const fix = (option == 1) || (option == 2)
        const include_45 = (option == 2)

        // Find the start and end index
        const start_index = find_start_index(MAG_Data[i].time)
        const end_index = find_end_index(MAG_Data[i].time)+1
        const num_samples = end_index - start_index

        // Get weighting based on bins
        const weights = get_weights(MAG_Data[i].bins.slice(start_index, end_index)).weights

        // Calculate average earth filed to match sensor to
        let ef_mean = { x:0.0, y:0.0, z:0.0 }
        for (let j = 0; j < num_samples; j++) {
            const data_index = start_index + j

            ef_mean.x += MAG_Data[i].orig.expected.x[data_index]
            ef_mean.y += MAG_Data[i].orig.expected.y[data_index]
            ef_mean.z += MAG_Data[i].orig.expected.z[data_index]
        }
        ef_mean.x /= num_samples
        ef_mean.y /= num_samples
        ef_mean.z /= num_samples

        let rotation = new Quaternion()

        // Try all rotations
        const last_rotation = 43
        let rot_error = []
        for (let rot = 0; rot <= last_rotation; rot++) {
            // Skip the weird ones
            if ((rot == 38) || (rot == 41)) {
                // ROTATION_ROLL_90_PITCH_68_yAW_293
                // ROTATION_PITCH_7
                continue
            }

            // Skip 45's if not enabled
            if (!include_45 && !right_angle_rotation(rot)) {
                continue
            }

            if (!rotation.from_rotation(rot)) {
                continue
            }

            // Rotate and take average
            let x = new Array(num_samples)
            let y = new Array(num_samples)
            let z = new Array(num_samples)
            let mean = { x:0.0, y:0.0, z:0.0 }
            for (let j = 0; j < num_samples; j++) {
                const data_index = start_index + j

                const tmp = rotation.rotate([MAG_Data[i].raw.x[data_index],
                                             MAG_Data[i].raw.y[data_index],
                                             MAG_Data[i].raw.z[data_index]])
    
                x[j] = tmp[0]
                y[j] = tmp[1]
                z[j] = tmp[2]

                mean.x += x[j]
                mean.y += y[j]
                mean.z += z[j]
            }
            mean.x /= num_samples
            mean.y /= num_samples
            mean.z /= num_samples

            const offsets = { 
                x: ef_mean.x - mean.x,
                y: ef_mean.y - mean.y,
                z: ef_mean.z - mean.z
            }

            let error_sum = 0
            for (let j = 0; j < num_samples; j++) {
                const data_index = start_index + j

                error_sum += ((x[j] - MAG_Data[i].orig.expected.x[data_index] + offsets.x)**2 +
                              (y[j] - MAG_Data[i].orig.expected.y[data_index] + offsets.y)**2 +
                              (z[j] - MAG_Data[i].orig.expected.z[data_index] + offsets.z)**2) * weights[j]
            }

            rot_error.push({ rotation: rot, error: error_sum / num_samples })

        }

        rot_error.sort((a, b) => a.error - b.error);

        const first = rot_error[0]
        const second = rot_error[1]

        const is_correct = (first.rotation == MAG_Data[i].params.orientation)
        const cost_ratio = second.error / first.error

        // best error must be half that of next best to be sure
        const check_valid = cost_ratio > 2

        const correct_txt = is_correct ? "correct" : "incorrect"
        let txt = "Mag " + (i+1) + " " + correct_txt + " orientation " + get_rotation_name(MAG_Data[i].params.orientation)
        if (!is_correct) {
            txt += ", best orientation: " + get_rotation_name(first.rotation)
        }
        txt += ", second best orientation: " + get_rotation_name(second.rotation)
        txt += ", cost ratio: " + (cost_ratio).toFixed(2)
        console.log(txt)

        // Ordinal rotation
        MAG_Data[i].rotation = MAG_Data[i].params.orientation

        if (check_valid && !is_correct) {
            // Found incorrect rotation
            if (fix) {
                MAG_Data[i].rotation = first.rotation

            } else {
                // Warn user but do not fix
                alert(
                    "Mag " + (i+1) + " possible incorrect orientation: " + get_rotation_name(MAG_Data[i].params.orientation) + "\n" +
                    "Should be: " + get_rotation_name(first.rotation) + " ?\n" +
                    "Cost ratio: " + (cost_ratio).toFixed(2)
                )

            }

        }

        // Apply rotation
        let rot = new Quaternion()
        rot.from_rotation(MAG_Data[i].rotation)

        const len = MAG_Data[i].raw.x.length
        MAG_Data[i].rotated = { x: new Array(len), y: new Array(len), z: new Array(len) }
        for (let j = 0; j < len; j++) {
            const tmp = rot.rotate([ MAG_Data[i].raw.x[j],
                                     MAG_Data[i].raw.y[j],
                                     MAG_Data[i].raw.z[j] ])

            MAG_Data[i].rotated.x[j] = tmp[0]
            MAG_Data[i].rotated.y[j] = tmp[1]
            MAG_Data[i].rotated.z[j] = tmp[2]
        }

    }

    const end = performance.now();
    console.log(`Orientation check took: ${end - start} ms`);
}

function get_body_frame_ef(tilt_correction, yaw) {

    const len = tilt_correction.length

    ret = { x: new Array(len), y: new Array(len), z: new Array(len) }

    let q = new Quaternion()
    for (i = 0; i < len; i++) {

        q.from_euler(0, 0, -yaw[i])

        tmp = q.rotate(earth_field.vector)
        tmp = tilt_correction[i].inverted().rotate(tmp)

        ret.x[i] = tmp[0]
        ret.y[i] = tmp[1]
        ret.z[i] = tmp[2]

    }

    return ret
}

let source
function select_body_frame_attitude() {

    if (source != null) {
        // No need to re-calc
        return
    }

    for (const ef of body_frame_earth_field) {
        if (ef.select.checked) {
            source = ef
        }
    }
    if (source == null) {
        alert("No attitude source selected")
        throw new Error()
    }

    // Interpolate expected to logged compass and calculate error
    for (let i = 0; i < 3; i++) {
        if (MAG_Data[i] == null) {
            continue
        }

        // Spherical interpolation between arrays of quaternions
        function array_slerp(values, index, query_index) {

            const len = query_index.length
            let ret = { q1: new Array(len), q2: new Array(len), q3: new Array(len), q4: new Array(len) }

            const last_value_index = index.length - 1
            let interpolate_index = 0
            for (let i = 0; i < len; i++) {

                if (query_index[i] <= index[0]) {
                    // Before start
                    ret.q1[i] = values.q1[0]
                    ret.q2[i] = values.q2[0]
                    ret.q3[i] = values.q3[0]
                    ret.q4[i] = values.q4[0]

                    continue
                }
                if (query_index[i] >= index[last_value_index]) {
                    // After end
                    ret.q1[i] = values.q1[last_value_index]
                    ret.q2[i] = values.q2[last_value_index]
                    ret.q3[i] = values.q3[last_value_index]
                    ret.q4[i] = values.q4[last_value_index]
                    continue
                }

                // increment index until there is a point after the target
                for (interpolate_index; interpolate_index < last_value_index; interpolate_index++) {
                    if (query_index[i] < index[interpolate_index+1]) {
                        const ratio = (query_index[i] - index[interpolate_index]) / (index[interpolate_index+1] - index[interpolate_index])

                        // Create A and C quaternions
                        const a = { q1: values.q1[interpolate_index],     q2: values.q2[interpolate_index],     q3: values.q3[interpolate_index],     q4: values.q4[interpolate_index]}
                        const c = { q1: values.q1[interpolate_index + 1], q2: values.q2[interpolate_index + 1], q3: values.q3[interpolate_index + 1], q4: values.q4[interpolate_index + 1]}

                        // interpolate
                        const b = slerp(a, c, ratio)

                        ret.q1[i] = b.q1
                        ret.q2[i] = b.q2
                        ret.q3[i] = b.q3
                        ret.q4[i] = b.q4

                        break
                    }
                }

            }
            return ret
        }

        MAG_Data[i].quaternion = array_slerp(source.quaternion, source.quaternion.time, MAG_Data[i].time)

        // Get yaw from quaternion for comparison later
        let quat = new Quaternion()
        let quat_yaw = new Quaternion()
        const len = MAG_Data[i].quaternion.q1.length
        MAG_Data[i].quaternion.yaw = new Array(len)
        MAG_Data[i].tilt_correction = new Array(len)
        for (let j = 0; j < len; j++) {

            // Populate quaternion
            quat.q1 = MAG_Data[i].quaternion.q1[j]
            quat.q2 = MAG_Data[i].quaternion.q2[j]
            quat.q3 = MAG_Data[i].quaternion.q3[j]
            quat.q4 = MAG_Data[i].quaternion.q4[j]

            // Record original yaw
            const yaw = quat.get_euler_yaw()
            MAG_Data[i].quaternion.yaw[j] = yaw

            // Quaternion to remove yaw
            quat_yaw.from_euler(0, 0, -yaw)

            // Remove yaw from attitude to provide a tilt correction rotation
            MAG_Data[i].tilt_correction[j] = quat_yaw.mul(quat)

        }

        // Rotate earth field into body frame
        MAG_Data[i].orig.expected = get_body_frame_ef(MAG_Data[i].tilt_correction, MAG_Data[i].quaternion.yaw)

        // Error between existing calibration and expected
        MAG_Data[i].orig.error = calc_error(MAG_Data[i].orig.expected, MAG_Data[i].orig)

        // Yaw estimate from existing calibration
        MAG_Data[i].orig.yaw = get_yaw(MAG_Data[i].orig, MAG_Data[i].tilt_correction)
    }
}

function fit() {

    const start = performance.now()

    // Run fit
    for (let i = 0; i < 3; i++) {
        if (MAG_Data[i] == null) {
            continue
        }

        // Find the start and end index
        const start_index = find_start_index(MAG_Data[i].time)
        const end_index = find_end_index(MAG_Data[i].time)+1
        const num_samples = end_index - start_index

        // Get weighting based on bins
        const weight_obj = get_weights(MAG_Data[i].bins.slice(start_index, end_index))
        const weights = weight_obj.weights
        const sqrt_weight = array_sqrt(weights)

        // Tilt correction quaternions
        const tilt_correction = MAG_Data[i].tilt_correction.slice(start_index, end_index)

        // Update coverage graphic
        MAG_Data[i].coverage.value = weight_obj.coverage

        // Calculate original fit error for selected samples only
        let error_sum = 0
        for (let j = 0; j < num_samples; j++) {
            error_sum += weights[j] * MAG_Data[i].orig.error[start_index + j]**2 
        }
        MAG_Data[i].orig.mean_error = Math.sqrt(error_sum / num_samples)

        let rot = {}
        let orientation
        if (!MAG_Data[i].rotate) {
            // Use raw directly
            rot.x = MAG_Data[i].raw.x.slice(start_index, end_index)
            rot.y = MAG_Data[i].raw.y.slice(start_index, end_index)
            rot.z = MAG_Data[i].raw.z.slice(start_index, end_index)

            // Original orientation
            orientation = MAG_Data[i].params.orientation

        } else {
            // use rotation corrected
            rot.x = MAG_Data[i].rotated.x.slice(start_index, end_index)
            rot.y = MAG_Data[i].rotated.y.slice(start_index, end_index)
            rot.z = MAG_Data[i].rotated.z.slice(start_index, end_index)

            // New orientation (possibly)
            orientation = MAG_Data[i].rotation

        }

        // Solve in the form Ax = B
        let A = new mlMatrix.Matrix(num_samples*3, 12)
        let B = new mlMatrix.Matrix(num_samples*3, 1)

        function setup_iron(A, row, colum, x, y, z) {

            const x_row = row + 0
            const y_row = row + 1
            const z_row = row + 2

            // Diagonal 1
            A.data[x_row][colum] = x
            A.data[y_row][colum] = 0.0
            A.data[z_row][colum] = 0.0

            // Diagonal 2
            colum++
            A.data[x_row][colum] = 0.0
            A.data[y_row][colum] = y
            A.data[z_row][colum] = 0.0

            // Diagonal 3
            colum++
            A.data[x_row][colum] = 0.0
            A.data[y_row][colum] = 0.0
            A.data[z_row][colum] = z

            // Off Diagonal 1
            colum++
            A.data[x_row][colum] = y
            A.data[y_row][colum] = x
            A.data[z_row][colum] = 0.0

            // Off Diagonal 2
            colum++
            A.data[x_row][colum] = z
            A.data[y_row][colum] = 0.0
            A.data[z_row][colum] = x

            // Off Diagonal 3
            colum++
            A.data[x_row][colum] = 0.0
            A.data[y_row][colum] = z
            A.data[z_row][colum] = y
        }

        function setup_offsets(A, row, colum, weight) {

            const x_row = row + 0
            const y_row = row + 1
            const z_row = row + 2

            // Offset 1
            A.data[x_row][colum] = weight
            A.data[y_row][colum] = 0.0
            A.data[z_row][colum] = 0.0

            // Offset 2
            colum++
            A.data[x_row][colum] = 0.0
            A.data[y_row][colum] = weight
            A.data[z_row][colum] = 0.0

            // Offset 3
            colum++
            A.data[x_row][colum] = 0.0
            A.data[y_row][colum] = 0.0
            A.data[z_row][colum] = weight

        }

        function setup_motor(A, row, colum, val) {

            const x_row = row + 0
            const y_row = row + 1
            const z_row = row + 2

            // Motor 1
            A.data[x_row][colum] = val
            A.data[y_row][colum] = 0.0
            A.data[z_row][colum] = 0.0

            // Motor 2
            colum++
            A.data[x_row][colum] = 0.0
            A.data[y_row][colum] = val
            A.data[z_row][colum] = 0.0

            // Motor 3
            colum++
            A.data[x_row][colum] = 0.0
            A.data[y_row][colum] = 0.0
            A.data[z_row][colum] = val

        }

        function setup_scale(A, row, colum, x, y, z) {

            const x_row = row + 0
            const y_row = row + 1
            const z_row = row + 2

            // Scale
            A.data[x_row][colum] = x
            A.data[y_row][colum] = y
            A.data[z_row][colum] = z

        }


        for (let fit of MAG_Data[i].fits) {

            let mot_val
            if (fit.value != null) {
                mot_val = fit.value.slice(start_index, end_index)
            }

            function evaluate_fit(expected, params) {

                function params_valid(params) {
                    function check_range(val, range) {
                        return (val > range[0]) && (val < range[1])
                    }

                    let ret = true
                    for (let i = 0; i < 3; i++) {
                        ret &= check_range(params.offsets[i], offsets_range)
                        ret &= check_range(params.diagonals[i], diagonals_range)
                        ret &= check_range(params.off_diagonals[i], off_diagonals_range)
                    }
                    ret &= check_range(params.scale, scale_range)
                    return ret
                }



                // Populate any unset params with defaults
                if (params.fit_type == null) {
                    params.fit_type = 0
                }
                if (params.diagonals == null) {
                    params.diagonals = [1.0, 1.0, 1.0]
                }
                if (params.off_diagonals == null) {
                    params.off_diagonals = [0.0, 0.0, 0.0]
                }
                if (params.scale == null) {
                    params.scale = 1.0
                }
                if (params.motor == null) {
                    params.motor = [0.0, 0.0, 0.0]
                }
                params.orientation = orientation

                // Check param ranges
                let ret = { params, valid: params_valid(params) }

                ret.valid = true

                if (!ret.valid) {
                    return ret
                }

                apply_params(ret, rot, params, mot_val)
                ret.error = calc_error(expected, ret)

                // Calculate error for selected samples only
                let error_sum = 0
                for (let j = 0; j < num_samples; j++) {
                    error_sum += weights[j] * ret.error[j]**2
                }
                ret.mean_error = Math.sqrt(error_sum / num_samples)

                ret.yaw = get_yaw(ret, MAG_Data[i].tilt_correction)

                return ret
            }

            function iterate_solution(fit_inst, B_population_fun, correction_fun, name, extract_param_fun) {
                // Start with the estimate from the original log
                fit_inst.expected = {
                    x: MAG_Data[i].orig.expected.x.slice(start_index, end_index),
                    y: MAG_Data[i].orig.expected.y.slice(start_index, end_index),
                    z: MAG_Data[i].orig.expected.z.slice(start_index, end_index),
                }

                // Pre-decompose matrix A
                const QR_A = new mlMatrix.QrDecomposition(A)

                let converged = false
                let params
                let yaw

                // Run limited number of iterations
                let k
                for (k = 0; k < 1000; k++) {

                    // Setup B matrix
                    B_population_fun(fit_inst.expected)

                    // Solve
                    params = QR_A.solve(B)

                    // Get corrected earth field
                    const corrected = correction_fun(params)

                    // Update yaw
                    const new_yaw = get_yaw(corrected, tilt_correction)

                    // Check for convergence
                    if (yaw != null) {
                        let max_change = 0
                        for (let j = 0; j < num_samples; j++) {
                            max_change = Math.max(max_change, Math.abs(wrap_PI(yaw[j] - new_yaw[j])))
                        }
                        converged = max_change < 0.001 * (Math.PI / 180.0)
                    }
                    yaw = new_yaw

                    // Recalculate expected
                    fit_inst.expected = get_body_frame_ef(tilt_correction, yaw)

                    if (converged) {
                        // Were done
                        break
                    }

                }

                if (converged) {
                    console.log("Mag " + (i+1) + " " + name + " " + fit.name + " converged after " + (k+1))
                } else {
                    console.log("Mag " + (i+1) + " " + name + " " + fit.name + " failed to converged after " + k)
                }

                Object.assign(fit_inst, evaluate_fit(fit_inst.expected, extract_param_fun(params)))
                //fit_inst.valid &= converged

            }


            const fit_mot = fit.value != null

            // Just fitting offsets, possibly with motor correction
            A.columns = fit_mot ? 6 : 3

            // A matrix
            for (let j = 0; j < num_samples; j++) {
                const index = j*3
                setup_offsets(A, j*3, 0, sqrt_weight[j])

                if (fit_mot) {
                    setup_motor(A, j*3, 3, mot_val[j] * sqrt_weight[j])
                }
            }

            function populate_B_offsets(expected) {
                for (let j = 0; j < num_samples; j++) {
                    const index = j*3
                    B.data[index+0][0] = (expected.x[j] - rot.x[j]) * sqrt_weight[j]
                    B.data[index+1][0] = (expected.y[j] - rot.y[j]) * sqrt_weight[j]
                    B.data[index+2][0] = (expected.z[j] - rot.z[j]) * sqrt_weight[j]
                }
            }

            function correct_offsets(params) {
                let corrected = {
                    x: array_offset(rot.x, params.get(0,0)),
                    y: array_offset(rot.y, params.get(1,0)),
                    z: array_offset(rot.z, params.get(2,0))
                }
                if (fit_mot) {
                    corrected.x = array_add(corrected.x, array_scale(fit.value, params.get(3,0)))
                    corrected.y = array_add(corrected.y, array_scale(fit.value, params.get(4,0)))
                    corrected.z = array_add(corrected.z, array_scale(fit.value, params.get(5,0)))
                }
                return corrected
            }

            function extract_offsets(params) {
                // Extract params
                let offsets = [ params.get(0,0), params.get(1,0), params.get(2,0) ]

                let motor
                if (fit_mot) {
                    motor = [ params.get(3,0), params.get(4,0), params.get(5,0) ]
                }

                return { offsets, motor, fit_type: fit.type }
            }

            // Solve
            iterate_solution(fit.offsets, populate_B_offsets, correct_offsets, "offsets", extract_offsets)


            // Just fitting offsets and scale, possibly with motor correction
            A.columns = fit_mot ? 7 : 4

            // Offsets already in column 0,1,2
            // Add scale and motor
            for (let j = 0; j < num_samples; j++) {
                const index = j*3

                setup_scale(A, index, 3, rot.x[j] * sqrt_weight[j], rot.y[j] * sqrt_weight[j], rot.z[j] * sqrt_weight[j])

                if (fit_mot) {
                    setup_motor(A, index, 4, mot_val[j] * sqrt_weight[j])
                }
            }

            function populate_B(expected) {
                for (let j = 0; j < num_samples; j++) {
                    const index = j*3
                    B.data[index+0][0] = expected.x[j] * sqrt_weight[j]
                    B.data[index+1][0] = expected.y[j] * sqrt_weight[j]
                    B.data[index+2][0] = expected.z[j] * sqrt_weight[j]
                }
            }

            function correct(params) {
                // Apply calibration
                const corrected_mat = A.mmul(params)
                let corrected = {
                    x: new Array(num_samples),
                    y: new Array(num_samples),
                    z: new Array(num_samples)
                }
                for (let j = 0; j < num_samples; j++) {
                    const index = j*3
                    corrected.x[j] = corrected_mat.data[index+0][0] / sqrt_weight[j]
                    corrected.y[j] = corrected_mat.data[index+1][0] / sqrt_weight[j]
                    corrected.z[j] = corrected_mat.data[index+2][0] / sqrt_weight[j]
                }
                return corrected
            }

            function extract_offsets_and_scale(params) {
                let scale = params.get(3,0)

                // Remove scale from offsets
                let offsets = array_scale([ params.get(0,0), params.get(1,0), params.get(2,0) ], 1 / scale)

                let motor
                if (fit_mot) {
                    motor = [ params.get(4,0), params.get(5,0), params.get(6,0) ]
                }

                return { offsets, scale, motor, fit_type: fit.type }
            }

            iterate_solution(fit.scale, populate_B, correct, "scale", extract_offsets_and_scale)

            // Fitting offsets and iron matrix, possibly with motor correction

            // Adjust size of A matrix depending if full mot fit is being done
            A.columns = fit_mot ? 12 : 9

            for (let j = 0; j < num_samples; j++) {
                const index = j*3

                setup_iron(A, index, 3, rot.x[j] * sqrt_weight[j], rot.y[j] * sqrt_weight[j], rot.z[j] * sqrt_weight[j])

                if (fit_mot) {
                    setup_motor(A, index, 9, mot_val[j] * sqrt_weight[j])
                }

            }

            function extract_offsets_and_iron(params) {
                // Extract params
                let diagonals =     [ params.get(3,0), params.get(4,0), params.get(5,0) ]
                let off_diagonals = [ params.get(6,0), params.get(7,0), params.get(8,0) ]

                // Remove iron correction from offsets
                const iron = new mlMatrix.Matrix([
                    [diagonals[0],     off_diagonals[0], off_diagonals[1]],
                    [off_diagonals[0], diagonals[1],     off_diagonals[2]], 
                    [off_diagonals[1], off_diagonals[2], diagonals[2]]
                ])
                const uncorrected_offsets = new mlMatrix.Matrix([[ params.get(0,0), params.get(1,0), params.get(2,0) ]])
                let offsets = Array.from(uncorrected_offsets.mmul(mlMatrix.inverse(iron)).data[0])

                // Normalize iron matrix into scale param
                let scale = array_mean(diagonals)
                diagonals = array_scale(diagonals, 1 / scale)
                off_diagonals = array_scale(off_diagonals, 1 / scale)

                let motor
                if (fit_mot) {
                    motor = [ params.get(9,0), params.get(10,0), params.get(11,0) ]
                }

                return { offsets, scale, diagonals, off_diagonals, motor, fit_type: fit.type }
            }

            iterate_solution(fit.iron, populate_B, correct, "iron", extract_offsets_and_iron)


            // Disable selection of invalid fits
            // select the first valid none motor fit by default
            let show_fit
            for (const key of Object.keys(fit_types)) {
                fit[key].show.disabled = !fit[key].valid
                if (fit[key].valid && (show_fit == null)) {
                    show_fit = fit[key].show
                }
            }
            if (show_fit && (fit.type == 0)) {
                show_fit.checked = true
            }

        }

    }

    const end = performance.now();
    console.log(`Fit took: ${end - start} ms`);

}

function calc_error(A, B) {
    const len = A.x.length
    let ret = new Array(len)
    for (let i = 0; i < len; i++) {
        ret[i] = Math.sqrt((A.x[i] - B.x[i])**2 + (A.y[i] - B.y[i])**2 + (A.z[i] - B.z[i])**2)
    }
    return ret
}

function extractLatLon(log) {
  var Lat, Lng
  if (("ORGN" in log.messageTypes) && ("instances" in log.messageTypes.ORGN) && (0 in log.messageTypes.ORGN.instances)) {
    const ORGN = log.get_instance("ORGN", 0)
    Lat = ORGN.Lat[ORGN.Lat.length-1] * 10**-7
    Lng = ORGN.Lng[ORGN.Lng.length-1] * 10**-7
    return [Lat, Lng]
  }
  console.warn("no ORGN message found")
  if ("POS" in log.messageTypes) {
    const POS = log.get("POS")
    Lat = POS.Lat[POS.Lat.length-1] * 10**-7
    Lng = POS.Lng[POS.Lng.length-1] * 10**-7
    return [Lat, Lng]
  }
  return [Lat, Lng]
}

// Enable/disable calculate and save params button
function set_need_calc(b) {
    document.getElementById('calculate').disabled = !b
    document.getElementById('SaveParams').disabled = b
}

function add_attitude_source(quaternion, name) {

    // Add check box for this attitude source
    let section = document.getElementById("ATTITUDE")

    let radio = document.createElement("input")
    radio.setAttribute('type', 'radio')
    radio.setAttribute('id', "ATTITUDE" + name)
    radio.setAttribute('name', "attitude_source")
    radio.disabled = false

    // Clear selected source and enable re-calc
    radio.addEventListener('change', function() { 
        set_need_calc(true)
        source = null
    })

    let label = document.createElement("label")
    label.setAttribute('for', "ATTITUDE" + name)
    label.innerHTML = name

    section.appendChild(radio)
    section.appendChild(label)
    section.appendChild(document.createElement("br"))

    return { quaternion, name, select: radio }

}

var MAG_Data
var fits
var body_frame_earth_field
var earth_field
async function load(log_file) {

    // Make sure imports are fully loaded before starting
    // This is needed when called from "open in"
    await import_done

    let log = new DataflashParser()
    log.processData(log_file, [])

    if (!("MAG" in log.messageTypes) || !("instances" in log.messageTypes.MAG)) {
        alert("No compass data in log")
        return
    }

    // micro seconds to seconds helpers
    const US2S = 1 / 1000000
    function TimeUS_to_seconds(TimeUS) {
        return array_scale(TimeUS, US2S)
    }

    // Plot flight data from log
    if ("ATT" in log.messageTypes) {
        const ATT_time = TimeUS_to_seconds(log.get("ATT", "TimeUS"))
        flight_data.data[0].x = ATT_time
        flight_data.data[0].y = log.get("ATT", "Roll")

        flight_data.data[1].x = ATT_time
        flight_data.data[1].y = log.get("ATT", "Pitch")
    } else {
        flight_data.data[0].x = null
        flight_data.data[0].y = null
        flight_data.data[1].x = null
        flight_data.data[1].y = null
    }

    if ("RATE" in log.messageTypes) {
        flight_data.data[2].x = TimeUS_to_seconds(log.get("RATE", "TimeUS"))
        flight_data.data[2].y = log.get("RATE", "AOut")
    } else {
        flight_data.data[2].x = null
        flight_data.data[2].y = null
    }

    if ("POS" in log.messageTypes) {
        flight_data.data[3].x = TimeUS_to_seconds(log.get("POS", "TimeUS"))
        flight_data.data[3].y = log.get("POS", "RelHomeAlt")
    } else {
        flight_data.data[3].x = null
        flight_data.data[3].y = null
    }

    Plotly.redraw("FlightData")

    // html helper
    function half_gap() {
        let hr = document.createElement("hr")
        hr.style.visibility = "hidden"
        hr.style.margin = "5px"
        return hr
    }

    MAG_Data = []

    // Helper to create tool tip returning a element
    function add_tip(parent, text, html) {
        parent.appendChild(document.createTextNode(" "))

        let img = document.createElement("img")
        parent.appendChild(img)

        img.src = "../images/question-circle.svg"
        img.style.width = "1em"
        img.style.verticalAlign = "bottom"
        img.setAttribute('data-tippy-content', text)
        img.setAttribute('data-tippy-maxWidth', '750px')

        if (html === true) {
            img.setAttribute('data-tippy-allowHTML', 'true')
        }

        tippy(img)
        return img
    }

    const PARM = log.get("PARM")
    function get_param(name, allow_change) {
        return get_param_value(PARM, name, allow_change)
    }

    // Get MAG data
    MAG_Data.start_time = null
    MAG_Data.end_time = null
    for (let i = 0; i < 3; i++) {

        // Clear section
        let name = "MAG" + i
        let info = document.getElementById(name)
        info.replaceChildren()

        if (!(i in log.messageTypes.MAG.instances)) {
            info.appendChild(document.createTextNode("Not found"))
            document.getElementById("MAG_" + (i + 1) + "_PARAM_INFO").replaceChildren(document.createTextNode("Not found"))
            continue
        }

        const MAG_msg = log.get_instance("MAG", i)

        // Load data from log
        MAG_Data[i] = { orig: { x: Array.from(MAG_msg.MagX),
                                y: Array.from(MAG_msg.MagY),
                                z: Array.from(MAG_msg.MagZ)},
                        time: TimeUS_to_seconds(MAG_msg.TimeUS),
                        names: get_compass_param_names(i+1),
                        fits: [],
                        param_selection: [] }

        // Set start and end times
        MAG_Data[i].start_time = MAG_Data[i].time[0]
        MAG_Data[i].end_time = MAG_Data[i].time[MAG_Data[i].time.length - 1]

        MAG_Data.start_time = (MAG_Data.start_time == null) ? MAG_Data[i].start_time : Math.min(MAG_Data.start_time, MAG_Data[i].start_time)
        MAG_Data.end_time = (MAG_Data.end_time == null) ? MAG_Data[i].end_time : Math.max(MAG_Data.end_time, MAG_Data[i].end_time)

        // Get param values
        MAG_Data[i].params = { offsets:  [ get_param(MAG_Data[i].names.offsets[0]),
                                           get_param(MAG_Data[i].names.offsets[1]),
                                           get_param(MAG_Data[i].names.offsets[2])],
                               diagonals: [ get_param(MAG_Data[i].names.diagonals[0]),
                                            get_param(MAG_Data[i].names.diagonals[1]),
                                            get_param(MAG_Data[i].names.diagonals[2])],
                               off_diagonals: [ get_param(MAG_Data[i].names.off_diagonals[0]),
                                                get_param(MAG_Data[i].names.off_diagonals[1]),
                                                get_param(MAG_Data[i].names.off_diagonals[2])],
                               scale: get_param(MAG_Data[i].names.scale),
                               motor: [ get_param(MAG_Data[i].names.motor[0]),
                                        get_param(MAG_Data[i].names.motor[1]),
                                        get_param(MAG_Data[i].names.motor[2])],
                               id: get_param(MAG_Data[i].names.id),
                               use: get_param(MAG_Data[i].names.use),
                               external: get_param(MAG_Data[i].names.external),
                               orientation: get_param(MAG_Data[i].names.orientation) }


        // Print some device info, offset is first param in fieldset
        const id = decode_devid(MAG_Data[i].params.id, DEVICE_TYPE_COMPASS)
        if (id != null) {
            if (id.bus_type_index == 3) {
                // DroneCAN
                info.appendChild(document.createTextNode(id.bus_type + " bus: " + id.bus + " node id: " + id.address + " sensor: " + id.sensor_id))
            } else {
                info.appendChild(document.createTextNode(id.name + " via " + id.bus_type))
            }
        }

        info.appendChild(half_gap())

        info.appendChild(document.createTextNode("Use: " + (MAG_Data[i].params.use ? "\u2705" : "\u274C")))
        info.appendChild(document.createTextNode(", "))
        info.appendChild(document.createTextNode("External: " + ((MAG_Data[i].params.external > 0) ? "\u2705" : "\u274C")))
        info.appendChild(document.createTextNode(", "))
        info.appendChild(document.createTextNode("Health: " + (array_all_equal(MAG_msg.Health, 1) ? "\u2705" : "\u274C")))

        info.appendChild(half_gap())

        info.appendChild(document.createTextNode("Coverage: "))
        MAG_Data[i].coverage = document.createElement("progress")
        info.appendChild(MAG_Data[i].coverage)
        add_tip(info, "This represents how many vehicle orientations are present in the log, higher coverage give more confidence in the results. It is out of all possible orientations (including upside-down), greater than 30% is good.")

        info.appendChild(half_gap())

        // Remove calibration to get raw values

        // Subtract compass-motor compensation
        let x = array_sub(MAG_Data[i].orig.x, Array.from(MAG_msg.MOX))
        let y = array_sub(MAG_Data[i].orig.y, Array.from(MAG_msg.MOY))
        let z = array_sub(MAG_Data[i].orig.z, Array.from(MAG_msg.MOZ))

        // Remove iron correction
        if (!array_all_equal(MAG_Data[i].params.diagonals, 0.0)) {

            // Invert iron correction matrix
            const iron = new mlMatrix.Matrix([
                [ MAG_Data[i].params.diagonals[0],     MAG_Data[i].params.off_diagonals[0], MAG_Data[i].params.off_diagonals[1] ],
                [ MAG_Data[i].params.off_diagonals[0], MAG_Data[i].params.diagonals[1],     MAG_Data[i].params.off_diagonals[2] ],
                [ MAG_Data[i].params.off_diagonals[1], MAG_Data[i].params.off_diagonals[2], MAG_Data[i].params.diagonals[2] ]
            ])
            const inv_iron = mlMatrix.inverse(iron)

            // Vectorized multiplication
            const corrected_x = array_add(array_add( array_scale(x, inv_iron.get(0,0)), array_scale(y, inv_iron.get(0,1))), array_scale(z, inv_iron.get(0,2)) )
            const corrected_y = array_add(array_add( array_scale(x, inv_iron.get(1,0)), array_scale(y, inv_iron.get(1,1))), array_scale(z, inv_iron.get(1,2)) )
            const corrected_z = array_add(array_add( array_scale(x, inv_iron.get(2,0)), array_scale(y, inv_iron.get(2,1))), array_scale(z, inv_iron.get(2,2)) )

            x = corrected_x; y = corrected_y; z = corrected_z
        }

        // Remove scale factor, if valid
        if (scale_valid(MAG_Data[i].params.scale)) {
            const inv_scale = 1 / MAG_Data[i].params.scale
            x = array_scale(x, inv_scale)
            y = array_scale(y, inv_scale)
            z = array_scale(z, inv_scale)
        }

        // remove offsets
        x = array_sub(x, Array.from(MAG_msg.OfsX))
        y = array_sub(y, Array.from(MAG_msg.OfsY))
        z = array_sub(z, Array.from(MAG_msg.OfsZ))

        // Rotate external compasses back into raw sensor frame
        let rotation = new Quaternion()
        let rotate = false
        if ((MAG_Data[i].params.external != 0) && rotation.from_rotation(MAG_Data[i].params.orientation)) {
            rotation.invert()
            const len = x.length
            for (let j = 0; j < len; j++) {
                const tmp = rotation.rotate([x[j], y[j], z[j]])
    
                x[j] = tmp[0]
                y[j] = tmp[1]
                z[j] = tmp[2]
            }
            rotate = true
        }

        MAG_Data[i].raw = { x: x, y: y, z: z }
        MAG_Data[i].rotate = rotate
    }

    // Set start and end time
    document.getElementById("TimeStart").value = MAG_Data.start_time
    document.getElementById("TimeEnd").value = MAG_Data.end_time

    // Assume constant earth field
    // Use origin msg
    // Use last EKF origin for earth field
    var [Lat, Lng] = extractLatLon(log)
    earth_field = expected_earth_field_lat_lon(Lat, Lng)
    if (earth_field == null) {
        alert("Could not get earth field for Lat: " + Lat + " Lng: " + Lng)
        return
    }
    console.log("EF: " + earth_field.vector[0] + ", " + earth_field.vector[1] + ", " + earth_field.vector[2] + " at Lat: " + Lat + " Lng: " + Lng)

    // Workout which attitude source to use, Note that this is not clever enough to deal with primary changing in flight
    const EKF_TYPE = get_param("AHRS_EKF_TYPE")

    // Load various attitude sources and calculate body frame earth field

    // Clear attitude selection options
    let attitude_select = document.getElementById("ATTITUDE")
    attitude_select.replaceChildren(attitude_select.children[0])

    body_frame_earth_field = []
    source = null

    if ("AHR2" in log.messageTypes) {

        const quaternion = {
            time: TimeUS_to_seconds(log.get("AHR2", "TimeUS")),
            q1: Array.from(log.get("AHR2", "Q1")),
            q2: Array.from(log.get("AHR2", "Q2")),
            q3: Array.from(log.get("AHR2", "Q3")),
            q4: Array.from(log.get("AHR2", "Q4"))
        }

        let field = add_attitude_source(quaternion, "DCM")
        if (EKF_TYPE == 0) {
            field.select.checked = true
        }
 
        body_frame_earth_field.push(field)
    }

    if (("NKQ" in log.messageTypes) && ("instances" in log.messageTypes.NKQ) && (0 in log.messageTypes.NKQ.instances)) {

        const quaternion = {
            time: TimeUS_to_seconds(log.get_instance("NKQ", 0, "TimeUS")),
            q1: Array.from(log.get_instance("NKQ", 0, "Q1")),
            q2: Array.from(log.get_instance("NKQ", 0, "Q2")),
            q3: Array.from(log.get_instance("NKQ", 0, "Q3")),
            q4: Array.from(log.get_instance("NKQ", 0, "Q4"))
        }

        let field = add_attitude_source(quaternion, "EKF 2 IMU 1")
        if (EKF_TYPE == 2) {
            field.select.checked = true
        }

        body_frame_earth_field.push(field)
    }

    if (("XKQ" in log.messageTypes) && ("instances" in log.messageTypes.XKQ)) {

        var primary = 0
        const EKF3_PRIMARY = get_param("EK3_PRIMARY")
        if (EKF3_PRIMARY != null) {
            primary = EKF3_PRIMARY
        }

        if (primary in log.messageTypes.XKQ.instances) {

            const quaternion = { 
                time: TimeUS_to_seconds(log.get_instance("XKQ", primary, "TimeUS")),
                q1: Array.from(log.get_instance("XKQ", primary, "Q1")),
                q2: Array.from(log.get_instance("XKQ", primary, "Q2")),
                q3: Array.from(log.get_instance("XKQ", primary, "Q3")),
                q4: Array.from(log.get_instance("XKQ", primary, "Q4"))
            }

            let field = add_attitude_source(quaternion, "EKF 3 IMU " + (primary + 1))
            if (EKF_TYPE == 3) {
                field.select.checked = true
            }
            body_frame_earth_field.push(field)
        }
    }

    if (body_frame_earth_field.length == 0) {
        alert("Unknown attitude source")
        return
    } else if (body_frame_earth_field.length == 1) {
        // Only one item, select it and disable
        body_frame_earth_field[0].select.checked = true
        body_frame_earth_field[0].select.disabled = true
    }

    // Add interference sources
    motor_comp.data = []

    // No compass motor fit
    for (let i = 0; i < 3; i++) {
        if (MAG_Data[i] == null) {
            continue
        }

        const tip_text = "Calibrations with no motor compensation.<ul><li>Offsets: calibration of only X Y Z offset parameters.</li><li>Offsets and scale: calibration of X Y Z offsets and single scale factor.</li><li>Offsets and iron: calibration of X Y Z offsets and iron compensation matrix diagonals and off-diagonals.</li></ul> If a calibration cannot be selected the tool was unable to find a valid solution, a longer flight with better coverage will give a better chance of finding a solution."

        MAG_Data[i].fits.push({
            value: null,
            type: 0,
            name: "No motor comp",
            tip: { text: tip_text, html: true },
            offsets: {},
            scale: {},
            iron: {}
        })
    }

    if (("BAT" in log.messageTypes) && ("instances" in log.messageTypes.BAT)) {
        for (let i = 0; i < 1; i++) {
            if (!(i in log.messageTypes.BAT.instances)) {
                continue
            }
            const value = Array.from(log.get_instance("BAT", i, "Curr"))
            if (array_all_NaN(value) || array_all_equal(value, 0)) {
                // Battery does not support current
                continue
            }
            const time = TimeUS_to_seconds(log.get_instance("BAT", i, "TimeUS"))
            for (let j = 0; j < 3; j++) {
                if (MAG_Data[j] == null) {
                    continue
                }
                const name =  "Battery " + (i+1) + " current"
                MAG_Data[j].fits.push({
                    value: linear_interp(value, time, MAG_Data[i].time),
                    type: 2,
                    name: name,
                    tip: { text: "Calibrations with motor compensation from " + name + ". Ensure battery monitor is calibrated and functioning correctly."},
                    offsets: {},
                    scale: {},
                    iron: {}
                })
            }

            // Add to motor comp plot
            motor_comp.data.push({
                mode: "lines",
                name: (i + 1).toFixed(),
                meta: "Battery " + (i + 1).toFixed(),
                legendgroup: 1,
                legendgrouptitle: { text: "Battery current" },
                hovertemplate: "<extra></extra>%{meta}<br>%{x:.2f} s<br>%{y:.2f} A",
                x: time,
                y: value
            })
        }
    }

    // Redraw motor comp plot
    Plotly.newPlot("motor_comp", motor_comp.data, motor_comp.layout, {displaylogo: false});

    // Hide if there is no data
    const hide_motor_comp = motor_comp.data.length == 0
    let plot = document.getElementById("motor_comp")
    plot.hidden = hide_motor_comp
    plot.previousElementSibling.hidden = hide_motor_comp

    // Add button for each fit
    for (let i = 0; i < 3; i++) {
        if (MAG_Data[i] == null) {
            continue
        }
        let name = "MAG" + i
        let section = document.getElementById(name)

        // Fieldset to contain parameter change options
        let param_fieldset = document.createElement("fieldset")
        section.appendChild(param_fieldset)

        let param_legend = document.createElement("legend")
        param_legend.innerHTML = "Parameter changes"
        add_tip(param_legend, "Include parameter change in saved parameter file. Use sensor allows the use parameter to be changed. Orientation allows the tool to fix incorrect orientations to the nearest 90 deg or 45 deg orientations. The tool will alert you if the orientation may be incorrect.")

        param_fieldset.appendChild(param_legend)

        function setup_radio(parent, type, label_txt, value) {
            const id = name + type + label_txt
            let radio = document.createElement("input")
            radio.setAttribute('type', 'radio')
            radio.setAttribute('id', id)
            radio.setAttribute('name', name + type)
            radio.setAttribute('value', value)

            let label = document.createElement("label")
            label.setAttribute('for', id)
            label.innerHTML = label_txt

            parent.appendChild(radio)
            parent.appendChild(label)

            return radio
        }

        param_fieldset.appendChild(document.createTextNode("Use sensor: "))
        setup_radio(param_fieldset, "use", "No change", 0).checked = true
        param_fieldset.appendChild(document.createTextNode(", "))
        setup_radio(param_fieldset, "use", "Use", 1)
        param_fieldset.appendChild(document.createTextNode(", "))
        setup_radio(param_fieldset, "use", "Don't use", 2)

        param_fieldset.appendChild(half_gap())

        param_fieldset.appendChild(document.createTextNode("Orientation:  "))
        
        let orientation = setup_radio(param_fieldset, "orientation", "Check", 0)
        orientation.addEventListener('change', function() { loading_call(() => { calculate(); }) } )
        orientation.checked = true
        orientation.disabled = !MAG_Data[i].rotate

        param_fieldset.appendChild(document.createTextNode(", "))
        orientation = setup_radio(param_fieldset, "orientation", "Fix 90", 1)
        orientation.addEventListener('change', function() { loading_call(() => { calculate(); }) } )
        orientation.disabled = !MAG_Data[i].rotate

        param_fieldset.appendChild(document.createTextNode(", "))
        orientation = setup_radio(param_fieldset, "orientation", "Fix 45", 2)
        orientation.addEventListener('change', function() { loading_call(() => { calculate(); }) } )
        orientation.disabled = !MAG_Data[i].rotate

        // Fieldset to contain calibration selection options
        let cal_fieldset = document.createElement("fieldset")
        section.appendChild(cal_fieldset)

        let cal_legend = document.createElement("legend")
        cal_legend.innerHTML = "Calibrations"
        add_tip(cal_legend, 'Select calibrations to be shown on plots, the last calibration selected will be saved when "Save Parameters" is clicked.')
        cal_fieldset.appendChild(cal_legend)



        function setup_check(parent, type, fit) {
            const id = name + type + fit
            let check = document.createElement("input")
            check.setAttribute('type', 'checkbox')
            check.setAttribute('id', id)
            check.addEventListener('change', function() { update_hidden(this) } )

    
            let label = document.createElement("label")
            label.setAttribute('for', id)
            label.innerHTML = type

            parent.appendChild(check)
            parent.appendChild(label)
            parent.appendChild(document.createElement("br"))

            return check
        }
        MAG_Data[i].orig.show = setup_check(cal_fieldset, "Existing", "")
        MAG_Data[i].orig.show.checked = true
        MAG_Data[i].orig.show.style.margin = "3px 3px 9px 20px"

        for (let j = 0; j < MAG_Data[i].fits.length; j++) {
            let fieldset = document.createElement("fieldset")

            let legend = document.createElement("legend")
            legend.innerHTML = MAG_Data[i].fits[j].name
            add_tip(legend, MAG_Data[i].fits[j].tip.text, MAG_Data[i].fits[j].tip.html)

            fieldset.appendChild(legend)

            for (const [key, value] of Object.entries(fit_types)) {
                MAG_Data[i].fits[j][key].show = setup_check(fieldset, value, MAG_Data[i].fits[j].name)
            }

            cal_fieldset.appendChild(fieldset)
        }
    }

    calculate()

}

// Update flight data range and enable calculate when time range inputs are updated
function time_range_changed() {

    flight_data.layout.xaxis.range = [ parseFloat(document.getElementById("TimeStart").value),
                                       parseFloat(document.getElementById("TimeEnd").value)]
    flight_data.layout.xaxis.autorange = false
    Plotly.redraw("FlightData")

    set_need_calc(true)
}
