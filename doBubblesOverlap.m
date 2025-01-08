function [isoverlap,ovV,iov,jov]=doBubblesOverlap(BUBBLES)

iov=0;
jov=0;
positions=BUBBLES.XYZ;
bubble_diam=BUBBLES.DIAM;
nbubbles_sim=numel(bubble_diam);

isoverlap=false;
ovV=0;

overlapVolume = @(r1, r2, d) (pi * (r1 + r2 - d)^2 * (d^2 + 2*d*r1 - 3*r1^2 + 2*d*r2 - 3*r2^2 + 6*r1*r2)) / (12 * d);


for i = 1:nbubbles_sim


    for j = i+1:nbubbles_sim


        distance = sqrt(sum((positions(i,:) - positions(j,:)).^2));

        % Check if spheres overlap
        if distance < (bubble_diam(i)/2 + bubble_diam(j)/2)

            ovV=overlapVolume(bubble_diam(i)/2,bubble_diam(j)/2,distance)./((4/3)*pi*((bubble_diam(i)/2).^3 + (bubble_diam(j)/2).^3));

            isoverlap=true;

            iov=i;
            jov=j;

            break

        end
    end
end


end