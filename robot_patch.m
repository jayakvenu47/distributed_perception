function [ patch_data ] = robot_patch()
    % Make it facing 0 rads
    robot_width = 0.6;
    robot_height = 0.8; 
    wheel_width = 0.1; 
    wheel_height = 0.4; 
    led_size = 0.1; 
    
    % Helper functions to generate vertex coordinates for a centered
    % rectangle and a helper function to shift a rectangle.
    rectangle = @(w, h) [w/2 h/2 1; -w/2 h/2 1; -w/2 -h/2 1; w/2 -h/2 1];
    shift = @(r, x, y) r + repmat([x, y, 0], size(r, 1), 1);
    
    % Create vertices for body, wheel, and led.
    body = rectangle(robot_width, robot_height);
    wheel = rectangle(wheel_width, wheel_height);
    led = rectangle(led_size, led_size);
    
    % Use pre-generated vertices and shift them around to create a robot
    left_wheel = shift(wheel, -(robot_width + wheel_width)/2, 0);
    right_wheel = shift(wheel, (robot_width + wheel_width)/2, 0);
    main_led = shift(led,  0, robot_height/4);

    
    % Putting all the robot vertices together
    vertices = [
     body ; 
     left_wheel; 
     right_wheel;
     main_led
    ];

    % Only color the body of the robot.  Everything else is black.
    % blue 85 157 205 contrast color 205 132 85
    colors = [
     [85 157 205]/255; 
     [1 1 1] *0.3;
     [1 1 1] *0.3;
     1 1 1 
    ];

    % It tells the patch function which vertices to connect.
    faces = repmat([1 2 3 4 1], 4, 1);
    
    for i = 2:4
       faces(i, :) = faces(i, :) + (i-1)*4;
    end
    
   patch_data = []; 
   patch_data.vertices = vertices;
   patch_data.colors = colors;
   patch_data.faces = faces;
end