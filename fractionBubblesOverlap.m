function overlapfraction=fractionBubblesOverlap(positions,bubble_diam)

nbubbles_sim=numel(bubble_diam);

if nbubbles_sim==0
    overlapfraction=0;
else

noverlap=0;

for i = 1:nbubbles_sim


    for j = i+1:nbubbles_sim


        distance = sqrt(sum((positions(i,:) - positions(j,:)).^2));

        % Check if spheres overlap
        if distance < (bubble_diam(i)/2 + bubble_diam(j)/2)

            noverlap = noverlap + 1;

        end

    end

end

overlapfraction = noverlap./nbubbles_sim;

end

end