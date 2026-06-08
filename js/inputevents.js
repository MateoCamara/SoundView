var previous_move_x = 0;
var previous_move_y = 0;
var start_mouse_down_time = 0;
var play_on_hover = false;
var wheel_zoom = true;

function onWheel(event){
    // canvas.removeEventListener("mousemove", onMouseMove, false);
    var e = window.event || event; // old IE support
    var delta = e.deltaY;

    if (wheel_zoom) {
        zoom_factor -= delta * 0.0005; //delta_zoom / ((w + h)/2);
        if (zoom_factor < min_zoom) {
            zoom_factor = min_zoom;
        } else if (zoom_factor > max_zoom){
            zoom_factor = max_zoom;
        }
    }
}


function onKeyDown(event){
    if (event.altKey){
        // enable play on hover mode
        play_on_hover = true;
        canvas.addEventListener("mousemove", onMouseMove, false);
    }
}

function onKeyUp(event){
    if (!event.altKey){
        // disable play on hover mode
        play_on_hover = false;
        canvas.removeEventListener("mousemove", onMouseMove, false);
    }
}

function onMouseDown(event){
    var d = new Date();
    start_mouse_down_time = d.getTime();
    previous_move_x = event.pageX;
    previous_move_y = event.pageY;
    canvas.addEventListener("mousemove", onMouseMove, false);
}

function onMouseUp(event){
    canvas.removeEventListener("mousemove", onMouseMove, false);
    previous_move_x = 0;
    previous_move_y = 0;

    var d = new Date();
    current_click_time = d.getTime();
    if ((current_click_time - start_mouse_down_time) < 100){
        click_x = event.pageX;
        click_y = event.pageY;
        var norm_click_x = 0;
        var norm_click_y = 0;
        [norm_click_x, norm_click_y] = displayCoordsToNormCoords(click_x, click_y);
        checkSelectSound(norm_click_x, norm_click_y);
    }
}

function onMouseOut(event){
    canvas.removeEventListener("mousemove", onMouseMove, false);
    previous_move_x = 0;
    previous_move_y = 0;
}

function onMouseMove(event){
    var zoomMode = event.metaKey;
    var rotateMode = false; // Disable rotation as it will make axis meaningless event.shiftKey;
    // var playMode = play_on_hover; // event.metaKey;
    var playMode = false; // event.metaKey;

    // Compute translation delta
    delta_x = event.pageX - previous_move_x;
    delta_y = event.pageY - previous_move_y;
    previous_move_x = event.pageX;
    previous_move_y = event.pageY;
    
    if (zoomMode) {
        zoom_factor -= delta_y * 0.05; //delta_zoom / ((w + h)/2);
        if (zoom_factor < min_zoom) {
            zoom_factor = min_zoom;
        } else if (zoom_factor > max_zoom){
            zoom_factor = max_zoom;
        }
    }

    if (rotateMode){
        rotation_degrees += delta_y * 0.01; //delta_zoom / ((w + h)/2);
    }

    if ((!zoomMode) && (!rotateMode) && (!playMode)) {
        var sin = Math.sin(-rotation_degrees);
        var cos = Math.cos(-rotation_degrees);
        var delta_x_r = delta_x * cos - delta_y * sin;
        var delta_y_r = delta_x * sin + delta_y * cos;
        center_x -= delta_x_r / (disp_scale * zoom_factor);
        center_y -= delta_y_r / (disp_scale * zoom_factor);  
    }

    if (playMode){
        var pos_x = 0;
        var pos_y = 0;
        [pos_x, pos_y] = displayCoordsToNormCoords(event.pageX, event.pageY);
        // throttle(function() {
        //     checkSelectSound(pos_x, pos_y)
        // }, 250);
        checkSelectSound(pos_x, pos_y);
    }
    
    // Stop after X seconds of inactivity
    // cancel timeout 
    // Set time out
}

function throttle(fn, threshhold, scope) {
    threshhold || (threshhold = 250);
    var last,
        deferTimer;
    return function () {
        var context = scope || this;

        var now = +new Date,
            args = arguments;
        if (last && now < last + threshhold) {
            // hold on to it
            clearTimeout(deferTimer);
            deferTimer = setTimeout(function () {
                last = now;
                fn.apply(context, args);
            }, threshhold);
        } else {
            last = now;
            fn.apply(context, args);
        }
    };
}

