function BUBBLEOUT=combineMoveOverlappingBubbles(BUBBLES,maxcoord)

%bubblediam_thresh = 15e-6;
%nloops_thresh=inf;

% BUBBLES_small_index = BUBBLES.DIAM < bubblediam_thresh;
%
% bubbles_small_diam = BUBBLES.DIAM(BUBBLES_small_index);
% bubbles_small_xyz = BUBBLES.XYZ(BUBBLES_small_index,:);

% positions=BUBBLES.XYZ(~BUBBLES_small_index,:);
% bubble_diam=BUBBLES.DIAM(~BUBBLES_small_index);

positions=BUBBLES.XYZ;
bubble_diam=BUBBLES.DIAM;

bubble_vol=(4/3)*pi*(bubble_diam/2).^3;

nbubbles_sim=numel(bubble_diam);

checked=false(nbubbles_sim,nbubbles_sim);
overlapdetected=fractionBubblesOverlap(positions,bubble_diam)>0;

movefacbuffer=1e-3;
%maxsize_prob=5000e-6;
%minsize_prob=10e-6;
ntrysmovesphere_thresh=1e5;

x0 = 100e-6; % Midpoint
k = 1e6; % Steepness

sigmoid_probability = @(x) 1 ./ (1 + exp(-k * (x - x0)));

%distance_thresh=1e-6;

%zf=0.5e-6;
%mu=1;
%sigma=0.35; % N/m

%tdrain_thresh=1e-7;

%count_every_n_loops=10;

nloops=0;

%overlapfraction=1;


while overlapdetected

    overlapdetected=false;

    for i = 1:nbubbles_sim-1

        jvec=i+find(~checked(i,(i+1):nbubbles_sim));
        r1=bubble_diam(i)/2;


        for j = jvec

            checked(i,j)=true;

            distance = sqrt(sum((positions(i,:) - positions(j,:)).^2));
            r2=bubble_diam(j)/2;

            if distance < (r1 + r2)*(1-movefacbuffer)

                overlapdetected=true;

                %if distance<min([r1,r2])*movefacbuffer

                %   pcombine=1;

                %else

                %area_intersection = 2*pi*(1./(2*distance)).^2*(-distance+r1-r2)*(-distance-r1+r2)*(-distance+r1+r2)*(distance+r1+r2);
                %pcombine=(log10(area_intersection.^2 / (1/r1+1/r2))+26)/7;
                %pcombine=log10(2*min([r1,r2])/minsize_prob)/log10(maxsize_prob/minsize_prob);
                pcombine=sigmoid_probability(2*min([r1,r2]));

                if pcombine<0
                    pcombine=0;
                elseif pcombine>1
                    pcombine=1;
                end

                %end


                if pcombine>0 && pcombine<1

                    combinebubbles = binornd(1,pcombine);

                else

                    combinebubbles = logical(pcombine);

                end




                if ~combinebubbles
                    % Try to move bubble outside of other sphere along
                    % vector defined by their centers:
                    centerdistanceV = positions(i,:) - positions(j,:);
                    overlapV = (1+movefacbuffer).*centerdistanceV * ((r1+r2)-distance)./distance;
                    if bubble_diam(j)<bubble_diam(i)
                        kk=j;
                        shiftvec = -overlapV;
                    else
                        kk=i;
                        shiftvec = overlapV;
                    end

                    positionkknew = positions(kk,:) + shiftvec;

                    isnewsphereoutsideallotherspheres=~sum(sqrt(sum((positionkknew-positions([1:(kk-1),(kk+1):end],:)).^2,2))<(bubble_diam([1:(kk-1),(kk+1):end])/2+bubble_diam(kk)/2)*(1-movefacbuffer));

                    if isnewsphereoutsideallotherspheres
                        % if the sphere is outside all others, assign it to the new position,
                        % clear the checked variable
                        positions(kk,:)=positionkknew;
                        checked(kk,(kk+1):end)=false;
                        checked(:,kk)=false;
                        %disp('Found shifted')
                    else

                        % Else try to move the smaller sphere to other
                        % locations around the larger sphere

                        ntrysmovesphere=0;

                        while ~isnewsphereoutsideallotherspheres

                            if ntrysmovesphere>ntrysmovesphere_thresh
                                break
                            end
                            % Find new position for this sphere
                            phi = 2 * pi * rand;           % Azimuthal angle [0, 2*pi]
                            theta = acos(1-2*rand);        % Polar angle [0, pi] (uniform on the sphere)
                            x_ = (r1+r2) * (1+movefacbuffer) * sin(theta) * cos(phi);
                            y_ = (r1+r2) * (1+movefacbuffer) * sin(theta) * sin(phi);
                            z_ = (r1+r2) * (1+movefacbuffer) * cos(theta);
                            if x_>0 && x_<=maxcoord && y_>0 && y_<=maxcoord && z_>0 && z_<=maxcoord
                                positionkknew = positions(kk,:) + [x_, y_, z_]; % Translate the point by the reference point
                                isnewsphereoutsideallotherspheres=~sum(sqrt(sum((positionkknew-positions([1:(kk-1),(kk+1):end],:)).^2,2))<(bubble_diam([1:(kk-1),(kk+1):end])/2+bubble_diam(kk)/2)*(1-movefacbuffer));
                                ntrysmovesphere=ntrysmovesphere+1;
                            end

                        end
                    end

                    % Check if that worked
                    if isnewsphereoutsideallotherspheres
                        % if the sphere is outside all others, assign it to the new position,
                        % clear the checked variable
                        positions(kk,:)=positionkknew;
                        checked(kk,(kk+1):end)=false;
                        checked(:,kk)=false;
                        %disp('Found sphere')
                    else
                        % Else try to move the smaller sphere to other
                        % locations all over the 3d space


                        ntrysmovesphere=0;
                        while ~isnewsphereoutsideallotherspheres

                            if ntrysmovesphere>ntrysmovesphere_thresh
                                break
                            end
                            rand3_=rand(3,1);
                            x_ = rand3_(1)*maxcoord;
                            y_ = rand3_(2)*maxcoord;
                            z_ = rand3_(3)*maxcoord;
                            % Translate the point by the reference point
                            positionkknew = [x_, y_, z_];
                            %positionkknewall=[positionkknewall;positionkknew]
                            isnewsphereoutsideallotherspheres=~sum(sqrt(sum((positionkknew-positions([1:(kk-1),(kk+1):end],:)).^2,2))<(bubble_diam([1:(kk-1),(kk+1):end])/2+bubble_diam(kk)/2)*(1-movefacbuffer));
                            ntrysmovesphere=ntrysmovesphere+1;

                        end

                    end

                    % Checked if that worked
                    if isnewsphereoutsideallotherspheres
                        % if the sphere is outside all others, assign it to the new position,
                        % clear the checked variable
                        positions(kk,:)=positionkknew;
                        checked(kk,(kk+1):end)=false;
                        checked(:,kk)=false;
                        %disp('Found random')
                    else
                        % if this didn't work, combine it.
                        combinebubbles=true;
                        %disp('NOT FOUND')
                    end


                end

                if combinebubbles

                    new_vol=bubble_vol(i) + bubble_vol(j);
                    newpos=(positions(i,:) * bubble_vol(i) + positions(j,:) * bubble_vol(j)) / new_vol;

                    bubble_vol(i)=new_vol;
                    positions(i,:)=newpos;
                    bubble_diam(i)=2*(new_vol*(3/(4*pi)))^(1/3);
                    checked(i,(i+1):end)=false;
                    checked(:,i)=false;

                    bubble_vol(j)=[];
                    positions(j,:)=[];
                    bubble_diam(j)=[];
                    checked(j,:)=[];
                    checked(:,j)=[];

                    nbubbles_sim=nbubbles_sim-1;

                end

            end

            if overlapdetected
                break
            end

        end

        if overlapdetected
            break
        end


        % if nloops>=count_every_n_loops
        %     overlapfraction=fractionBubblesOverlap(positions,bubble_diam);
        %     nloops=0;
        % end

    end

    nloops=nloops+1;

end


%nloops
BUBBLEOUT=struct('XYZ',positions,'DIAM',bubble_diam);


% if numel(positions)==0
%     BUBBLEOUT=struct('XYZ',bubbles_small_xyz,'DIAM',bubbles_small_diam);
% else
%     BUBBLEOUT=struct('XYZ',cat(1,bubbles_small_xyz,positions),'DIAM',[bubbles_small_diam;bubble_diam]);
% end


end