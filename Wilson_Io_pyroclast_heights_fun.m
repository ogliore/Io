function Wilson_Io_pyroclast_heights_fun(total_SO2_mass_percentage_start)

%clearvars
%close all

docombinebubbles=true;

%pool=parpool;
%numWorkers = pool.NumWorkers;

doplots=false;
domovie=false;
dodebugging=false;

ndepthsteps=1e2;
bubble_sim_volume_m3 = 1e-12; %1e-12;
grid_size_npx=512; %512;
primary_bubble_size_microns=10;
primary_bubble_size_m=primary_bubble_size_microns*1e-6;

% combineOverlappingBubbles threshholds:
bubble_volfactor_thresh=1.5;
probability_to_combine_bubbles=0.1;
%bubblediam_thresh=primary_bubble_size_m*bubble_volfactor_thresh^(1/3);
%overlapfraction_thresh=0.50;

magmagasmixture_mass=1;
rho_t=2800; % kg/m^3

LL=load('RayleighPlesset_10um.mat');

total_SO2_mass_percentage_start=10;
maxdepth_meters=19e3;


bubble_sim_mass_initial = rho_t.*bubble_sim_volume_m3;

bubble_sim_mass_factor = bubble_sim_mass_initial./magmagasmixture_mass;

gas_vol_percentage_critical=75; % fragmentation

grid_volume=true(grid_size_npx,grid_size_npx,grid_size_npx);

bubble_sim_width = (bubble_sim_volume_m3).^(1/3);

bubble_sim_width_original=bubble_sim_width;

pixels_per_m=(grid_size_npx./bubble_sim_width);

pixels_per_m_original=pixels_per_m;

microns_per_pixel=1e6./pixels_per_m;

if domovie

    v = VideoWriter('bubblemovie.avi'); %, 'MPEG-4');
    open(v);  % Open the video file for writing

end

Tliq = 1700;
Tsol = 1420;
Tm=210;
sr = 1600;
L = 6e5;
E = 111.1;

Tb_from_n = @(n) -((n./100 - 1).*(Tliq.*sr + (L.*Tliq)./(Tliq - Tsol)) - (E.*Tm.*n)./100)./((E.*n)./100 - (n./100 - 1).*(sr + L./(Tliq - Tsol)));

if doplots
    nv__=linspace(0.1,30,1e4);
    plot(nv__,Tb_from_n(nv__),'LineWidth',3);
    xlabel('Mass fraction SO_2 in magma (%)','FontSize',20)
    ylabel('Bulk magma temperature (K)','FontSize',20)
    set(gca,'Linewidth',1.5,'FontSize',20)
    saveas(gca,'Magma_temp_vs_SO2.png')
    %keyboard
end

% Don't need the Iotherm for this calculation:
% A=importdata('Spencer2020_FIg2b_Iotherm.csv');
% Iotherm_temperature=A(:,1);
% Iotherm_r=A(:,2);
% Iotherm_depth=max(Iotherm_r)-Iotherm_r;
% [Iotherm_depth,sind]=sort(Iotherm_depth);
% Iotherm_temperature=Iotherm_temperature(sind);
% [Iotherm_depth,IC,~]=unique(Iotherm_depth);
% Iotherm_temperature=Iotherm_temperature(IC);

% Extrapolation value for high end = highest point?
%Io_temperature_from_depth = @(depth) interp1(Iotherm_depth*1e3,Iotherm_temperature,depth,'pchip');

g_Io=1.796; % m/s2

m_SO2=64.06e-3 ; % molecular mass of SO2 (kg/mol)

R=8.314462618;

rho_Io_crust=2792;


Io_pressure_from_depth = @(depth) rho_Io_crust .* g_Io .* depth;
Io_depth_from_pressure = @(pressure) pressure./(rho_Io_crust .* g_Io);

% Mysen and Popp (1980)
pressure_SO2_sol=[0,500,1e3,1.5e3,2e3,3e3]*1e6;
SO2_sol=[0.0,0.3,0.525,0.67,0.72,0.81];

% Li and Ripley (2009)
XH2O=0;
XFeO=0.136;
XTiO2=0.0098;
XCaO=0.122;
XSiO2=0.732;
T=1700;
Xs_from_T_P = @(T,P) 1e2*exp(-1.76-0.474*1e4/T-0.021*(P*1e-8)+5.559*XFeO+2.565*XTiO2+2.709*XCaO-3.192*XSiO2-3.049*XH2O);

Lionel_quartic = @(P) 6.265596*1e-4*(P*1e-6) - 1.294010e-8*(P*1e-6)^2 - 1.094052e-10*(P*1e-6)^3 + 2.469919e-14*(P*1e-6)^4;


if doplots
    ppv=logspace(3,8,1e3);
    %semilogx(ppv/1e6,Xs_from_T_P_allP(1700,ppv),pressure_SO2_sol./1e6,SO2_sol,'LineWidth',2)
    semilogx(ppv/1e6,Xs_from_T_P_allP(1700,ppv),'LineWidth',3)
    xlabel('Pressure (MPa)','FontSize',16)
    ylabel('S solubility (%)','FontSize',16)
    %legend('Li and Ripley (2009)','Mysen and Popp (1980)','Location','Best','FontSize',16)
    set(gca,'Linewidth',1.5,'FontSize',17.5)
    %xticks(logspace(0,5,5));
    %xticklabels(string(logspace(0,5,5)));
    saveas(gca,'SO2_solubility_vs_P_comparison_combined.png')
end


dv=linspace(5e3,0,1e6);
depthvalues=[linspace(maxdepth_meters,5e3,ndepthsteps-round(0.9*ndepthsteps)),interp1(Xs_from_T_P_allP(T,Io_pressure_from_depth(dv)),dv,linspace(Xs_from_T_P_allP(T,Io_pressure_from_depth(dv(1))),...
    Xs_from_T_P_allP(T,Io_pressure_from_depth(dv(end))),round(0.9*ndepthsteps)))];
dv=[];


SO2_sol_from_pressure = @(p_) interp1(pressure_SO2_sol,SO2_sol,p_);

magma_SO2_mass_percentage=nan(ndepthsteps,1);

magma_SO2_mass_percentage(1) = SO2_sol_from_pressure(Io_pressure_from_depth(depthvalues(1)));
exsolved_mass_percentage_start=total_SO2_mass_percentage_start-magma_SO2_mass_percentage(1);
gas_mass=(total_SO2_mass_percentage_start/100)*magmagasmixture_mass;
magma_mass=(1-total_SO2_mass_percentage_start/100).*magmagasmixture_mass;
magma_volume=magma_mass./rho_t;

gas_vol_percentage=nan(ndepthsteps,1);
gas_vol_percentage_grid=nan(ndepthsteps,1);
exsolved_mass_percentage=nan(ndepthsteps,1);
P=nan(ndepthsteps,1);
exsolved_SO2_mass=nan(ndepthsteps,1);
bubble_diameter_microns=nan(ndepthsteps,1);
exsolved_SO2_volume=nan(ndepthsteps,1);
ncombinedbubbles=nan(ndepthsteps,1);

exsolved_mass_percentage(1)=exsolved_mass_percentage_start;

P(1)=Io_pressure_from_depth(depthvalues(1));
T=Tb_from_n(total_SO2_mass_percentage_start);


max_bubble_factor=10;

max_bubble_diameter_microns=max_bubble_factor*primary_bubble_size_microns;

exsolved_SO2_mass(1)=magmagasmixture_mass*exsolved_mass_percentage(1)/100;
exsolved_SO2_moles_start=exsolved_SO2_mass(1)./m_SO2;
exsolved_SO2_volume_start=exsolved_SO2_moles_start*R*T(1)/P(1);
gas_vol_percentage(1)=100*exsolved_SO2_volume_start./(exsolved_SO2_volume_start+magma_volume);
bubble_diameter_microns_start=primary_bubble_size_microns;
bubble_diameter_microns(1)=bubble_diameter_microns_start;
exsolved_SO2_volume(1)=exsolved_SO2_volume_start;

bubble_step_size_microns=primary_bubble_size_microns*(1-mean((Io_pressure_from_depth(depthvalues(2:end))./Io_pressure_from_depth(depthvalues(1:end-1))).^(1./3)));

hist_smaller_fac=5;

volume_of_one_primary_bubble_um3=(4/3)*pi*(primary_bubble_size_microns/2)^3;
volume_of_one_primary_bubble_m3=volume_of_one_primary_bubble_um3*(1e-6)^3;
primary_bubble_size_nm=primary_bubble_size_microns*1e3;
histogram_step_size=bubble_step_size_microns./hist_smaller_fac;
bubble_diameters_edges= 0:histogram_step_size:max_bubble_diameter_microns;
bubble_diameters_mid=(bubble_diameters_edges(1:end-1)+(bubble_step_size_microns./hist_smaller_fac))';
nbubble_diameters_=numel(bubble_diameters_mid);
bubble_histogram=zeros(nbubble_diameters_,1);

[~, ~, binind_primary_bubble_size] = histcounts(primary_bubble_size_microns, bubble_diameters_edges);



BUBBLES=struct('XYZ',[],'DIAM',[]);
nbubbles=zeros(ndepthsteps,1);
new_exsolved_SO2_primary_bubble_number_in_sim_volume=zeros(ndepthsteps,1);

P=Io_pressure_from_depth(depthvalues);
new_solubility_total=Xs_from_T_P_allP(T,P);


for jj=1:2
    if jj==2
        %keyboard
        exsolved_SO2_mass(ii_break)=0;
        exsolved_mass_percentage(ii_break)=0;
        gas_vol_percentage(ii_break)=0;
        exsolved_SO2_volume(ii_break)=0;
        bubble_diameter_microns(ii_break)=bubble_diameter_microns_start;
        startind=ii_break+1;
    else
        startind=2;
        elapsedTime1=0; elapsedTime2=0;
    end
    for ii=startind:ndepthsteps

        new_exsolved_SO2_mass_percentage=magma_SO2_mass_percentage(ii-1)-new_solubility_total(ii);

        exsolved_mass_percentage(ii)=exsolved_mass_percentage(ii-1)+new_exsolved_SO2_mass_percentage;
        new_exsolved_SO2_mass=magmagasmixture_mass*new_exsolved_SO2_mass_percentage/100;
        exsolved_SO2_mass(ii)=exsolved_SO2_mass(ii-1)+new_exsolved_SO2_mass;
        new_exsolved_SO2_moles=new_exsolved_SO2_mass./m_SO2;
        new_exsolved_SO2_volume=new_exsolved_SO2_moles*R*T/P(ii);
        exsolved_SO2_volume(ii)=exsolved_SO2_volume(ii-1)*P(ii-1)/P(ii)+new_exsolved_SO2_volume;
        bubble_diameter_microns(ii)=bubble_diameter_microns_start.*(exsolved_SO2_volume(ii)./exsolved_SO2_volume_start).^(1./3);
        gas_vol_percentage(ii)=100*exsolved_SO2_volume(ii)./(exsolved_SO2_volume(ii)+magma_volume);
        magma_SO2_mass_percentage(ii)=new_solubility_total(ii);
        gas_vol_percentage_cut=gas_vol_percentage(ii);

        gas_vol_percentage_grid(ii)=0;

        if jj==2
            nbubbles(ii) = numel(BUBBLES.DIAM);

            if nbubbles(ii)>0
                % Increase size of each previous bubble
                % Rayleigh-Plesset Equation!
                % Pvisc = 4*eta*rbdot/rb % rb=BUBBLES.DIAM/2, rbdot = drb/dt
                % Psurften = 2*gamma_surface_tension/rb %
                % gamma_surface_tension = surface tension at magma-gas interface
                % Pbubble=Io_pressure_from_depth(depthvalues(ii))+Pvisc+Psurften;
                bubble_increase_diam_fac=interp1(LL.Pv,LL.RtRo,P(ii))./interp1(LL.Pv,LL.RtRo,P(ii-1));
                %bubble_increase_diam_fac=((T/P(ii))./(T/P(ii-1)))^(1./3);
                BUBBLES.DIAM = BUBBLES.DIAM.*bubble_increase_diam_fac;                
                %fprintf('Combining overlapping bubbles (volume increase): ')
                maxcoord=grid_size_npx./pixels_per_m;
                tic
                if docombinebubbles
                    BUBBLES=combineMoveOverlappingBubbles(BUBBLES,maxcoord);
                end
                elapsedTime1=toc;
                nbubbles(ii) = numel(BUBBLES.DIAM);
                %fprintf('%.1f s\n',elapsedTime1);

            end

            %if numel(BUBBLES.DIAM)>0 || (new_exsolved_SO2_mass_percentage>0 && exsolved_mass_percentage(ii)>0)



            if new_exsolved_SO2_mass_percentage>0 % even if exsolved so2 percentage is negative
                % Add new 10-micron bubbles:
                bubble_sim_volume_m3=bubble_sim_volume_m3.*(magma_volume+exsolved_SO2_volume(ii))./(magma_volume+exsolved_SO2_volume(ii-1));
                bubble_sim_width = (bubble_sim_volume_m3).^(1/3);
                pixels_per_m_last=pixels_per_m;
                pixels_per_m=(grid_size_npx./bubble_sim_width);
                microns_per_pixel=1e6./pixels_per_m;
                new_exsolved_SO2_primary_bubble_number=new_exsolved_SO2_volume./volume_of_one_primary_bubble_m3;
                %new_exsolved_SO2_primary_bubble_number_in_sim_volume=round(new_exsolved_SO2_primary_bubble_number*bubble_sim_volume_m3/(magma_volume+new_exsolved_SO2_volume+exsolved_SO2_volume(ii)));
                %new_exsolved_SO2_primary_bubble_number_in_sim_volume(ii)=round(new_exsolved_SO2_primary_bubble_number*bubble_sim_volume_m3/(magma_volume+exsolved_SO2_volume(ii)));
                new_exsolved_SO2_primary_bubble_number_in_sim_volume(ii)=round(new_exsolved_SO2_primary_bubble_number*bubble_sim_mass_factor);
                %rr = randperm(grid_size_npx, new_exsolved_SO2_primary_bubble_number_in_sim_volume);
                maxcoord=grid_size_npx./pixels_per_m;
                fgridvolume=find(grid_volume);
                fgridvolumer=randsample(fgridvolume,new_exsolved_SO2_primary_bubble_number_in_sim_volume(ii));
                [xyz1,xyz2,xyz3] = ind2sub(size(grid_volume),fgridvolumer);
                vv=linspace(1,512,10);
                subplot(2,3,1); histogram(xyz1,vv); subplot(2,3,2); histogram(xyz2,vv); subplot(2,3,3); histogram(xyz3,vv);
                vv2=vv./pixels_per_m;
                if numel(BUBBLES.DIAM)>0
                subplot(2,3,4); histogram(BUBBLES.XYZ(:,1),vv2); subplot(2,3,5); histogram(BUBBLES.XYZ(:,2),vv2); subplot(2,3,6); histogram(BUBBLES.XYZ(:,3),vv2);
                end
                %pause
                %xyz = cat(1,r',c',v');
                
                %OLD RANDOM SEED:
                %xyz = randi([1, grid_size_npx],3,new_exsolved_SO2_primary_bubble_number_in_sim_volume(ii));
                %xyz = confineRandomSeedsToMagma(grid_size_npx,grid_volume,xyz);


                % Exsolve bubbles only where no bubbles were before
                %[x, y, z] = fast_ind2sub(grid_size_npx, grid_size_npx, grid_size_npx, find(grid_volume))
                %[xmagma, ymagma, zmagma] = ind2sub(size(grid_volume), find(grid_volume));

                % rr = randperm(grid_size_npx, new_exsolved_SO2_primary_bubble_number_in_sim_volume(ii));
                % NEWPOSx = xmagma(rr);
                % NEWPOSy = ymagma(rr);
                % NEWPOSz = zmagma(rr);
                %BUBBLES.XYZ = [BUBBLES.XYZ ; [xyz(1,:)',xyz(2,:)',xyz(3,:)']./pixels_per_m];
                BUBBLES.XYZ  = [BUBBLES.XYZ.*pixels_per_m_last/pixels_per_m ;[xyz1,xyz2,xyz3]./pixels_per_m];
                %BUBBLES.XYZ = [BUBBLES.XYZ ; bubble_sim_width.*rand(new_exsolved_SO2_primary_bubble_number_in_sim_volume(ii),3)];
                BUBBLES.DIAM = [BUBBLES.DIAM.*pixels_per_m_last/pixels_per_m; repmat(primary_bubble_size_m,new_exsolved_SO2_primary_bubble_number_in_sim_volume(ii),1)];
                %fprintf('Combining overlapping bubbles (new bubbles): ')
                tic
                if dodebugging
                    fprintf('%.2e -> ',double(sum((4/3)*pi*(BUBBLES.DIAM(:)/2).^3)))
                end
                if docombinebubbles
                    BUBBLES=combineMoveOverlappingBubbles(BUBBLES,maxcoord);
                end
                if dodebugging
                    fprintf('%.2e, ',  double(sum((4/3)*pi*(BUBBLES.DIAM(:)/2).^3)))
                    [Btf,ovV,iov,jov]=doBubblesOverlap(BUBBLES);
                    if Btf
                        fprintf('OVERLAP: (%d,%d): %.3f%%\n',iov,jov,100*ovV);
                    else
                        fprintf('no overlap\n');
                    end
                end
                elapsedTime2 = toc;
                nbubbles(ii)=numel(BUBBLES.DIAM);
                %fprintf('%.1f s\n',elapsedTime2);

            end


            if nbubbles(ii)>0
                grid_volume=true(grid_size_npx,grid_size_npx,grid_size_npx);
                for mm=1:nbubbles(ii)
                    % Calculate the bounding box of the sphere (limit to within matrix bounds)
                    x_min = max(1, floor(pixels_per_m.*(BUBBLES.XYZ(mm,1) - 0.5*BUBBLES.DIAM(mm))));
                    x_max = min(grid_size_npx, ceil(pixels_per_m.*(BUBBLES.XYZ(mm,1) + 0.5*BUBBLES.DIAM(mm))));
                    y_min = max(1, floor(pixels_per_m.*(BUBBLES.XYZ(mm,2) - 0.5*BUBBLES.DIAM(mm))));
                    y_max = min(grid_size_npx, ceil(pixels_per_m.*(BUBBLES.XYZ(mm,2) + 0.5*BUBBLES.DIAM(mm))));
                    z_min = max(1, floor(pixels_per_m.*(BUBBLES.XYZ(mm,3) - 0.5*BUBBLES.DIAM(mm))));
                    z_max = min(grid_size_npx, ceil(pixels_per_m.*(BUBBLES.XYZ(mm,3) + 0.5*BUBBLES.DIAM(mm))));
                    [xGrid, yGrid, zGrid] = ndgrid(x_min:x_max, y_min:y_max, z_min:z_max);

                    distance_squared = (xGrid - pixels_per_m.*BUBBLES.XYZ(mm,1)).^2 + (yGrid - pixels_per_m.*BUBBLES.XYZ(mm,2)).^2 + (zGrid - pixels_per_m.*BUBBLES.XYZ(mm,3)).^2;
                    outside_bubble_is_true = distance_squared > (pixels_per_m.*0.5*BUBBLES.DIAM(mm)).^2;
                    %slice(double(outside_bubble_is_true), round(size(outside_bubble_is_true,1)/2), round(size(outside_bubble_is_true,2)/2), round(size(outside_bubble_is_true,3)/2)); axis equal; colorbar
                    grid_volume(x_min:x_max, y_min:y_max, z_min:z_max) = outside_bubble_is_true & grid_volume(x_min:x_max, y_min:y_max, z_min:z_max);
                    grid_size_new=floor(grid_size_npx*pixels_per_m./pixels_per_m_original);
                    grid_volume2=grid_volume(1:grid_size_new,1:grid_size_new,1:grid_size_new);

                end

            end

            gas_vol_percentage_grid(ii)=100*(1-sum(grid_volume(:))/grid_size_npx^3);

            if domovie & nbubbles(ii)>0
                f=figure('visible','off');
                p = patch(isosurface(grid_volume, 0.5));  % Extract surface at the boundary where the value is 0.5
                p.FaceColor = 'red';
                p.EdgeColor = 'none';
                daspect([1 1 1]);
                view(3);
                camlight; lighting gouraud;  % Add lighting for better visualization
                axis off
                % Adjust the axes to reduce whitespace
                ax = gca;  % Get the current axes
                outerpos = ax.OuterPosition;
                ti = ax.TightInset;
                left = outerpos(1) + ti(1);
                bottom = outerpos(2) + ti(2);
                ax_width = outerpos(3) - ti(1) - ti(3);
                ax_height = outerpos(4) - ti(2) - ti(4);
                ax.Position = [left, bottom, ax_width, ax_height];  % Set tight position for axes
                frame = getframe(gcf);
                writeVideo(v, frame);
            end

        end

        fprintf('D = %e (%d/%d), Sol = %.4e, Exs SO2 = %.4f%%, P = %.2e kPa, Gas = %.2f%%/%.2f%%, (%.1f/%.1f) s, %d Bs: %d new, max = %.2f um\n',...
            depthvalues(ii),ii,ndepthsteps,new_solubility_total(ii),exsolved_mass_percentage(ii),P(ii)./1e3,gas_vol_percentage_grid(ii),gas_vol_percentage(ii),elapsedTime1,elapsedTime2,numel(BUBBLES.DIAM),new_exsolved_SO2_primary_bubble_number_in_sim_volume(ii),max(BUBBLES.DIAM.*1e6));

        if (jj==1 && gas_vol_percentage(ii)>gas_vol_percentage_critical) || (jj==2 && gas_vol_percentage_grid(ii)>gas_vol_percentage_critical) 
            ii_break=ii;
            fprintf('n/mass%% = %.2f, Df1 = %.2f km, Pf1 = %.2f MPa, F1=%.4f\n',total_SO2_mass_percentage_start,depthvalues(ii_break)/1e3,P(ii_break)/1e6,(P(1)/P(ii_break)).^(1/3))
            break
        end

        histogram(BUBBLES.DIAM*1e6)
    end
end

if domovie
    close(v);
end

datetimenow=string(datetime('now', 'Format', 'dd-MM-yyyy-HH-mm-ss'));
datetimenow=datetimenow{1};



%%


% Label connected components of voids
rsz=2^3;
gvs= imresize3(grid_volume, size(grid_volume)/rsz,'cubic');
CC = bwconncomp(gvs, 26); % 26-connectivity for 3D
volumes = cellfun(@numel, CC.PixelIdxList);
diameters_um=1e6*2*(volumes.*(3/(4*pi))).^(1/3)/pixels_per_m;
histogram(diameters_um*rsz)

%%
npointsforhull=12;
nnshortest=returnClosestBubblePositions(BUBBLES.XYZ,BUBBLES.DIAM,npointsforhull);
interstitial_diameter_microns=nan(numel(BUBBLES.DIAM),1);


for ii=1:numel(BUBBLES.DIAM)    
    [~, convexhullV] = convhull(BUBBLES.XYZ(nnshortest(ii,:),1), BUBBLES.XYZ(nnshortest(ii,:),2), BUBBLES.XYZ(nnshortest(ii,:),3));
    sphereV=(4/3) * pi * (BUBBLES.DIAM(ii)/2).^3;
    neighboringspheresV = (4/3) * pi * (BUBBLES.DIAM(nnshortest(ii,:))/2).^3;
    interstitialV=convexhullV-sphereV-sum(0.5*neighboringspheresV);
    if interstitialV>0
        interstitial_diameter_microns(ii) = 1e6*2*(interstitialV*(3/(4*pi)))^(1/3);
    end
end

if doplot
    histogram(interstitial_diameter_microns)
    xlabel('Pyroclast diameter (microns)','FontSize',16);
    ylabel('Number','FontSize',16);
    set(gca,'Linewidth',1.5,'FontSize',17.5);
    saveas(gca,'Pyroclast_diameter_histogram.png')
end

vs = volshow(grid_volume);
obj = ancestor(vs,'figure','toplevel');
I = getframe(obj);
imwrite(I.cdata,'grid_volume.png');


interstitial_diameter_microns_5_95=prctile(interstitial_diameter_microns,[5,95]);

%%

save(sprintf('Wilson_Io_pyroclast_heights_%d_%s.mat',total_SO2_mass_percentage_start,datetimenow),'BUBBLES','grid_volume',...
    'total_SO2_mass_percentage_start','bubble_sim_volume_m3','depthvalues',...
    'bubble_sim_width','pixels_per_m','microns_per_pixel',...
    'volumes','diameters_um','CC',...
    'new_solubility_total','exsolved_mass_percentage','P','gas_vol_percentage_grid','gas_vol_percentage',...
    'new_exsolved_SO2_primary_bubble_number_in_sim_volume','interstitial_diameter_microns_5_95');



%%

%originalSize=size(grid_volume);

% Define the downsampling factor (block size)
% downsampleFactor = 2;
%
% newSizeX = originalSize(1) / downsampleFactor;
% newSizeY = originalSize(2) / downsampleFactor;
% newSizeZ = originalSize(3) / downsampleFactor;


% Reshape the matrix into blocks and apply the OR operation
%grid_volume_ds = squeeze( ...
%    any(reshape(~grid_volume, downsampleFactor, newSizeX, downsampleFactor, newSizeY, downsampleFactor, newSizeZ), [1 3 5]));

%ngrid_volume=~grid_volume;

% CC = bwconncomp(grid_volume,6);
%
% stats = regionprops3(CC, 'Volume', 'BoundingBox', 'Centroid', 'VoxelIdxList');
%
% % Visualize the regions
% figure;
% hold on;
%
% nregions=height(stats);
% ctbl = distinguishable_colors(nregions);
%
% % Loop through each region and visualize the 3D volume
% for i = 2:nregions
%     % Extract the indices of the voxels belonging to the current region
%     voxelIndices = stats.VoxelIdxList{i};
%
%     % Create a blank volume to hold the shape of the current region
%     regionVolume = false(size(~grid_volume));
%     regionVolume(voxelIndices) = true;
%
%     % Plot the 3D surface of the region
%     p = patch(isosurface(regionVolume, 0.5));
%     p.FaceColor = ctbl(i,:);  % Set color of the region surface
%     p.EdgeColor = 'none';  % Remove edges for a smoother look
% end
%
% % Format the 3D view
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');
% view(3);  % 3D view
% axis tight;
% camlight;
% lighting gouraud;  % Add lighting for better visualization
% hold off;
%
% interstitial_volumes=stats.Volume;
%
% interstitial_volumes(interstitial_volumes==max(interstitial_volumes))=[];
%
% figure
% interstitial_diameters_microns=microns_per_pixel.*2*(stats.Volume(2:end)*(3/(4*pi))).^(1/3);
% histogram(interstitial_diameters_microns,1:50)
% figure
% histogram(BUBBLES.DIAM*1e6)
%
%
%
% save('Wilson_Io_pyroclast_heights.mat')
%
%
