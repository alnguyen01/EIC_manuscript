% function [ alt_wat,alt_org_exp,alt_org_lay,alt_min,alt_ice_simple,alt_ice_simple_ONLY_ORG,alt_ice_simple_ONLY_ORG_MAX,defo,porosity_volume,alt_ice] = alt_from_deform_V3(deform,satfrac,Sandfrac,kRoot,mass_org,root_depth, org_depth,rho_om_max, poros_om, ABoVE_ALT)
function [ alt_wat,alt_org_exp,alt_org_lay,alt_min,defo] = alt_from_deform_V3(deform,satfrac,Sandfrac,kRoot,mass_org,root_depth, org_depth,rho_om_max, poros_om, ABoVE_ALT)
%UNTITLED Summary of this function goes here
% comment out function 1 for uncertainty in ALT -> defo; vice versa for
% uncertainty in deform -> ALT
%   Detailed explanation goes here

%!
%!=%======================================================================
%    SUBROUTINE alt_from_deform_V2(deform,satfrac, sandfrac,kroot,mass_org,&
%    root_depth, org_depth,rho_om_max, poros_om,alt_wat,alt_org_exp,alt_org_lay,alt_min)
%!=======================================================================
%! Calculate the active layer thickness given surface deformation assuming 
%%! different vertical profiles of porosity with depth
%! Version 2 is effective 11/20/14
%!
%! Modifications:
%%!  Kevin Schaefer made subroutine (8/30/10)
%!  Kevin Schaefer added organic layer (1/1/12)
%!  Kevin Schaefer added saturation fraction scaling (11/19/14)
% %!-----------------------------------------------------------------------
% %!
%     IMPLICIT NONE
% !
% ! Input variables
%     real deform     ! (cm) surface deformation
%     real satfrac    ! (-) saturation fraction of soil porosity
%     real sandfrac   ! (%) sand fraction of soil texture
%     real kRoot      ! (1/m) Exp const for decrease root-zone organic matter w/ depth
%     real mass_org   ! (kg/m2) mass of organic matter in top 1 m of soil
%     real root_depth ! (m) maximum rooting depth
%     real org_depth  ! (m) thickness of organic layer
%     real rho_om_max ! (kg/m3) maximum organic matter density
%     real poros_om   ! (-) soil porosity for pure organic soil
% !
% ! Output variables
%     real alt_wat     ! (cm) active layer thickness pure water
%     real alt_org_exp ! (cm) active layer thickness exponential mix of organic and mineral soil
%     real alt_org_lay ! (cm) active layer thickness organic soil layer
%     real alt_min     ! (cm) active layer thickness pure mineral soil
% !
% ! Local variables
%     integer iter    ! (-) iteration index
%     integer n_iter  ! (-) number iterations
%     real dz         ! (m) incremental active layer thickness
%     real depth      ! (m) current depth
%     real def_test   ! (cm) current deformation
%     real del_def    ! (cm) change in deformation
%     real def_max    ! (cm) deformation for layer of organic soil
%     real rho_water  ! (kg/m3) density of water
%     real rho_ice    ! (kg/m3) density of ice
%     real poros      ! (-) porosity for organic mineral soil mix
%     real poros_min  ! (-) mineral soil porosity
%     real orgfrac    ! (-) organic soil fraction
%     real fsat       ! (-) soil air space saturation scaling factor
%     real sat_max    ! (-) soil saturation fraction where fsat is a maximum
%     real sat_zero   ! (-) soil saturation fraction where fsat is zero
%     real sat_eff    ! (-) effective soil saturation fraction accounting for expansion into pore air
% !
% ! execution control variables
 %   logical :: flg_fsat=.true. ! (-) apply saturation fraction scaling factor
%!
%! assign ice and water densities
    rho_water=1000.;
    rho_ice=917.;
    size1=size(deform);
%!
%! vertical integration delta soil depth
%! tests indicate that this value produces minimal error
    dz=0.01;
%!
%! calculate mineral soil porosity
   poros_min=0.489-0.00126*Sandfrac;
%!
%! soil saturation scaling factor to account for air space
%! transfer saturation fraction to local variable 
%! to allow data assimilation of satfrac without stomping on it
%! bottom stop saturation since retrieval does not work for low saturation
   %if(flg_fsat) 
%!
%! scale saturation using Degesse et al. [2010] curve fit for 11.1% clay
%! curve fit has artifact peak at sat_max, which is less than 1.0
     sat_max=0.9629210   ;    
     if(satfrac>=sat_max) 
       fsat=1.;
     else
       fsat=-19.63107*satfrac.^2.+37.8063409*satfrac-17.2022565;
%!
%! bottom stop fsat to prevent unrealistically large ALT near sat_zero
       fsat=max(fsat, 0.1);
     end
%!
%! set lower limit on saturation based on zero deformation from Degesse et al. [2010] curve
     sat_zero=0.7372225;
     sat_eff=zeros(size(satfrac));
     %for i=1:size1(1)
         %for j=1:size1(2)
            sat_eff=max(satfrac,sat_zero);
         %end
     %end
%!
% scale saturation to account for expansion into pore air space
     sat_eff=sat_eff.*fsat;
  % else
%!
%! set lower limit of 0.2 on saturation (arbitrarily chosen)
     sat_zero=0.2;
     sat_eff=max(satfrac,sat_zero);
   %end
%!
%! make sure saturation does not exceed 1.0
    sat_eff=min(sat_eff,1.0);
    size1=size(deform);
%!-----------------------------------------------------------------------
%! alt_wat: ALT for pure water column
%!-----------------------------------------------------------------------
    %alt_wat=zeros(size1(1),size1(2));
    %alt_min=zeros(size1(1),size1(2));
    %alt_org_lay=zeros(size1(1),size1(2));
    alt_wat=rho_ice/(rho_water-rho_ice).*deform;
%!
%!-----------------------------------------------------------------------
%! alt_min: ALT for pure mineral soil
%!-----------------------------------------------------------------------
    alt_min=deform./sat_eff./poros_min.*rho_ice./(rho_water-rho_ice);
%!
%!-----------------------------------------------------------------------
%! alt_org_exp: exponential decrease in organic matter
%!-----------------------------------------------------------------------
%! numerical integration to calculate ALT for mix of organic and mineral soil 
    size1=size(deform);
    depth=zeros(size1(1),size1(2));
    def_test=zeros(size1(1),size1(2));
    del_def=zeros(size1(1),size1(2));
    def_max=zeros(size1(1),size1(2));
    dz=0.01 * ones(size1(1),size1(2));
   % n_iter=0;
    
    
    %Roger addon defs
    %root_depth=RootDepth;
    %rho_om_max=MaxOrganicMatter;
    %mass_org=OrganicMass;
    %poros_om=OrganicSoilPorosity;
    
  
        
     
%!
%! count total number iterations
    %deform(isnan(deform))=0;
    dz=0.01 * ones(size1(1),size1(2));
    n_iter=0;
    for i=1:size1(1)
        for j=1:size1(2)
            n_iter=0;
            while (abs(depth(i,j)-deform(i,j))>.001) && (n_iter<1000) %&& (isnan(deform(i,j))==1)
%! organic soil fraction
                orgfrac=kRoot*mass_org/(1.-exp(-kRoot*root_depth))*exp(-kRoot*(depth(i,j)+0.5*dz(i,j)))/rho_om_max;
                if orgfrac>1. 
                    orgfrac=1.;
                end
                if orgfrac<0. 
                    orgfrac=0.;
                end
                if depth(i,j)>root_depth 
                    orgfrac=0.;
                end
%!
%! soil porosity
                poros=(1.-orgfrac)*poros_min+orgfrac*poros_om;
                if poros>1. 
                    poros=1.;
                end
%!
%! change in deformation
                del_def(i,j)=poros*sat_eff(i,j)*(rho_water-rho_ice)/rho_ice*dz(i,j)*100.;
%!
%! check on current deformation and depth
                if(def_test(i,j)+del_def(i,j)>=deform(i,j)) %then ! deformation too big
                    dz(i,j)=0.5*dz(i,j);
                else %! continue with integration
                    def_test(i,j)=def_test(i,j)+del_def(i,j);
                    depth(i,j)=depth(i,j)+dz(i,j);
                end
%!
%! convergence test (deformation within 0.0001 cm)
%! tests indicate this is the best compromise between 
%! accuracy of estimated ALT and the number of iterations
                %if abs(def_test(i,j)-deform(i,j))>.001 
                n_iter=n_iter+1;
                %end
            end
         end
     end
    %end
%!
%! save ALT for mixed soil
    alt_org_exp=depth*100.;
%!%
%!-----------------------------------------------------------------------
%5! alt_org_lay: organic soil layer
%!-----------------------------------------------------------------------
%! organic soil fraction
    orgfrac=mass_org/org_depth/rho_om_max;
    orgfrac=max(orgfrac,0.);
    orgfrac=min(orgfrac,1.);
%!
%! soil porosity
    poros=(1.-orgfrac)*poros_min+orgfrac*poros_om;
    if poros>1. 
        poros=1.;
    end
%5!
alt_org_lay=zeros(size1(1),size1(2));
for i=1:size1(1)
    for j=1:size1(2)
%! max deformation due to organic layer
        def_max(i,j)=org_depth*sat_eff(i,j)*poros*(rho_water-rho_ice)/rho_ice*100.;
        if deform(i,j)<=def_max(i,j) 
            alt_org_lay(i,j)=deform(i,j)/sat_eff(i,j)/poros*rho_ice/(rho_water-rho_ice);
        else
            alt_org_lay(i,j)=def_max(i,j)./sat_eff(i,j)./poros.*rho_ice./(rho_water-rho_ice);
            alt_org_lay(i,j)=alt_org_lay(i,j)+(deform(i,j)-def_max(i,j))./sat_eff(i,j)./poros_min.*rho_ice./(rho_water-rho_ice);
        end
    end
end

    %RETURN
  

%end %this ends the function
%!-----------------------------------------------------------------------
%% alt_ground_ice_simple: ground ice layer beneath soil organic layer
%!-----------------------------------------------------------------------
%! assign ice porosity
   poros_ice = 0.0343; %Zhang, Qian, et al. 2022 ~for roughly same density as defined above
   poros_ice = 0.1275; %Zhang, Qian, et al. 2022 ~for at closest temp to 0 deg C (-5) provided 
%! organic soil fraction
    orgfrac=mass_org/org_depth/rho_om_max;
    orgfrac=max(orgfrac,0.);
    orgfrac=min(orgfrac,1.);
%!
%! soil porosity
    % poros=(1.-orgfrac)*poros_min+orgfrac*poros_om; %ORIGINAL!
    poros=(1.-orgfrac)*poros_ice+orgfrac*poros_om; %replacing poros_min with poros_ice
    if poros>1. 
        poros=1.;
    end
%5!
alt_ice_simple=zeros(size1(1),size1(2));
for i=1:size1(1)
    for j=1:size1(2)
%! max deformation due to organic layer
        def_max(i,j)=org_depth*sat_eff(i,j)*poros*(rho_water-rho_ice)/rho_ice*100.;
        if deform(i,j)<=def_max(i,j) 
            alt_ice_simple(i,j)=deform(i,j)/sat_eff(i,j)/poros*rho_ice/(rho_water-rho_ice);
            %NEXT COUPLE OF LINES ARE FOR VOLUMETRIC ICE CONTENT TO LIQUID
            %WATER!!!
            alt_ice_simple_ONLY_ORG(i,j)=deform(i,j)/sat_eff(i,j)/poros*rho_ice/(rho_water-rho_ice); %this is to be subtracted from the final result to show which pixels have org+ice
%only thing changed is poros_min -> poros_ice compared to org_layer model
        else
            alt_ice_simple(i,j)=def_max(i,j)./sat_eff(i,j)./poros.*rho_ice./(rho_water-rho_ice);
            alt_ice_simple_ONLY_ORG_MAX(i,j)=def_max(i,j)./sat_eff(i,j)./poros.*rho_ice./(rho_water-rho_ice); %this is to be subtracted from the result only org+ice to just get at ALT due to ice/pixel
            alt_ice_simple(i,j)=alt_ice_simple(i,j)+(deform(i,j)-def_max(i,j))./sat_eff(i,j)./poros_ice.*rho_ice./(rho_water-rho_ice);
        end
    end
end

    %RETURN

%% testing w/ ABoVE Map more...
% % Field Active Layer Thickness -> deformation
% %!-----------------------------------------------------------------------
% %! deform_exp: !!!!organic soil layer!!! -> deformation
% %inverse from eqn. 6 in Liu et al. 2012
% %!-----------------------------------------------------------------------
% dz = 1; %in cm
% 
% org_depth_cm = org_depth*100;
% org_depth_cm_index = org_depth_cm;
% Nmax = max(ABoVE_ALT, [],"all");
% 
% 
% defo=zeros(size1(1),size1(2));
% porosity_volume=zeros(size1(1),size1(2),floor(Nmax));
% 
% tic
% for i=1:size1(1)
%     for j=1:size1(2)
%         N = ABoVE_ALT(i,j)./dz; %number of steps
%         if isnan(N)
%             defo(i,j)=nan;
%         else
%             %set up porosity as a fxn of depth, step fxn
%             Porosity = zeros(1,floor(N)); %use fix to make sure N is an integer
%             for k = 1:org_depth_cm_index
%                 Porosity(1,k) = poros_om;
%             end
% 
%             for k = (org_depth_cm_index+1):floor(N)
%                 Porosity(1,k) = poros_min;
%             end
%             porosity_volume(i,j,1:floor(N)) = Porosity(1,1:floor(N));
%             Saturation = ones(1,floor(N)); %assume fully saturated soil
% 
%             % set up integration
%             defo0=0;
%             for m = 1:floor(N)
%                 defo0 = defo0 + Porosity(1,m) * Saturation(m) * dz * ((rho_water-rho_ice)/(rho_ice));
%             end
%             defo(i,j)=defo0;
%         end
%     end
% end
% toc

%% testing more w/ ABoVE Map but for exp organic instead of organic layer (see above)
% % Field Active Layer Thickness -> deformation
% %!-----------------------------------------------------------------------
% %! deform_exp: organic soil layer -> deformation
% %inverse from eqn. 6 in Liu et al. 2012
% %!-----------------------------------------------------------------------
%! alt_org_exp: exponential decrease in organic matter
%!-----------------------------------------------------------------------
%! numerical integration to calculate defo for ALT w/ mix of organic and mineral soil 
    Nmax = max(ABoVE_ALT, [],"all")./100; %meters
    size1=size(deform);


    defo=zeros(size1(1),size1(2));
    porosity_volume=zeros(size1(1),size1(2),floor(Nmax)); % was looking at this 

    dz=0.01; %meters, = 1 cm
    for i=1:size1(1)
        for j=1:size1(2)
            N = (ABoVE_ALT(i,j)./100)./dz; %number of steps, meters, = 1 cm
            if isnan(N)
                defo(i,j)=nan;
            else
                Porosity = zeros(1,floor(N));
                for k = 1:floor(N)
                    orgfrac=kRoot*mass_org/(1.-exp(-kRoot*root_depth))*exp(-kRoot*(k+0.5*dz))/rho_om_max;
                    if orgfrac>1. 
                        orgfrac=1.;
                    end
                    if orgfrac<0. 
                        orgfrac=0.;
                    end
                    if k>root_depth 
                        orgfrac=0.;
                    end
%!
%! soil porosity
                    Porosity(1,k)=(1.-orgfrac)*poros_min+orgfrac*poros_om;
                    if Porosity(1,k)>1. 
                        Porosity(1,k)=1.;
                    end
                end
                porosity_volume(i,j,1:floor(N)) = Porosity(1,1:floor(N));
                Saturation = ones(1,floor(N)); %assume fully saturated soil

                % set up integration
                defo0=0;
                for m = 1:floor(N)
                    defo0 = defo0 + Porosity(1,m) * Saturation(m) * dz * ((rho_water-rho_ice)/(rho_ice));
                end
                defo(i,j)=defo0.*100; %cm
            end
        end
    end
%!%
%% all ice (functionality is run the ABoVE map stuff to get EIC thickness
% then take the EIC thickness and run through here to get the effective
% "ALT" assuming certain things about porosity and saturation of an all
% ice layer

%! assign ice porosity
   poros_ice = 0.0343; %Zhang, Qian, et al. 2022 ~for roughly same density as defined above
   poros_ice = 0.1275; %Zhang, Qian, et al. 2022 ~for at closest temp to 0 deg C (-5) provided 

   %sat calculation
   vmc = 68.1721; %average VMC from all core measurements
   rho_s = 2.07;
   rho_dopn = 2.05;
   e = rho_s-rho_dopn/rho_dopn; %Dysli et al. 2003 scaling of supersaturation by a simple test
   rho_perm = 1.25; %somewhere between pure ice and no ice permafrost Kawasaki et al. 1983
   sat_ice = (vmc*rho_s)/(e*(rho_perm-vmc*rho_ice));

   alt_ice = deform./sat_ice./poros_ice.*rho_ice./(rho_water-rho_ice);


end