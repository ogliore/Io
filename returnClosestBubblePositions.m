function nnshortest=returnClosestBubblePositions(positions,bubble_diam,nn)

nbubbles_sim=numel(bubble_diam);

distance=nan(nbubbles_sim,nbubbles_sim);

for i = 1:nbubbles_sim


    for j = i+1:nbubbles_sim


        distance(i,j) = sqrt(sum((positions(i,:) - positions(j,:)).^2));


    end

end

nnshortest=nan(nbubbles_sim,nn);

for ii=1:nbubbles_sim
    [~,b]=sort(distance(ii,:));
    nnshortest(ii,:)=b(1:nn);
end