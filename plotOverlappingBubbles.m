
function plotOverlappingBubbles(bubble_diam,positions,i,j)

 if nargin==2
        colors = distinguishable_colors(size(positions, 1));
 end


figure;
hold on;
axis equal;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Spheres');
for mm = 1:size(positions, 1)

    % Generate sphere coordinates
    [x, y, z] = sphere(50); % Generate a sphere with finer resolution
    % Scale and shift the sphere
    x = bubble_diam(mm)/2 * x + positions(mm, 1);
    y = bubble_diam(mm)/2 * y + positions(mm, 2);
    z = bubble_diam(mm)/2 * z + positions(mm, 3);
    % Plot the sphere
    if nargin==4
        if mm==j
            surf(x, y, z, 'FaceColor', [1,0,0], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        elseif mm==i
            surf(x, y, z, 'FaceColor', [0,1,0], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        else
            surf(x, y, z, 'FaceColor', [0,0,1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        end
    end

    if nargin==2
            surf(x, y, z, 'FaceColor', colors(mm,:), 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        end
    end


end