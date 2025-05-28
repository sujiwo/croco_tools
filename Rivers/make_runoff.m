%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Create a CROCO run-off forcing file
%
%
%  Further Information:
%  http://www.crocoagrif.org
%
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%
%  July 2013: G. Cambon (IRD/LEGOS) & M. Herrmann (IRD/LEGOS)
%  July 2023 : G. Cambon (IRD/LOPS)
%
clear all
close all
%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%%
crocotools_param

%=========================================================================================
% Choose the grid level into which you ant to set up the runoffs

gridlevel=0;
if ( gridlevel == 0 ) 
    % #-> Parent / zoom #O 
    grdname  = [CROCO_files_dir,'croco_grd.nc'];
    rivname =  [CROCO_files_dir,'croco_runoff.nc'];
    clmname = [CROCO_files_dir,'croco_clm.nc'];    % <- climato file for runoff
else 
    % # -> Child / zoom #XX
    grdname  = [CROCO_files_dir,'croco_grd.nc.',num2str(gridlevel)];
    rivname =  [CROCO_files_dir,'croco_runoff.nc.',num2str(gridlevel)];
    clmname = [CROCO_files_dir,'croco_clm.nc.',num2str(gridlevel)]; % <- climato file for runoff                             
end

% Choose the monthly runoff forcing time and cycle in days
clim_run=0;

if (clim_run == 1)
    qbar_time=[15:30:345]; 
    qbar_cycle=360; 
else
    qbar_time=[15.2188:30.4375:350.0313];
    qbar_cycle=365.25;
end

%     - times and cycles for runoff conditions:
%           - clim_run = 1 % climato forcing experiments with climato calendar
%                     qbar_time=[15:30:365];
%                     qbar_cycle=360;
%
%           - clim_run = 0 % interanual forcing experiments with real calendar
%                     qbar_time=[15.2188:30.4375:350.0313];
%                     qbar_cycle=365.25;

%=========================================================================================
% Choose if you process variable tracer concentration(temp, salt, NO3, ...)

psource_ncfile_ts=0;

if psource_ncfile_ts
    psource_ncfile_ts_auto=1 ;
    psource_ncfile_ts_manual=0;
end

%        - psource_ncfile_ts = 0 => Constant analytical runoff tracers concentration no processing
%                            It reads analytical values in croco.in
%                            or use default value defined in
%                            analytical.F
%
%        - pource_ncfile_ts = 1  => Variable runoff tracers
%                                    concentration  processing is activated.
%                                   It needs the climatology
%                                   file created with make_clim.m 
%
%                                    In this case, either choose:
%                                      - psource_ts_auto : auto definition using
%                                                          the nearest point in the climatlogy file
%
%                                       - psource_ts_manual : manually definition the
%                                                             variable tracer concentration

%=========================================================================================
% Fancy plots

plotting=1;
plotting_zoom=0;
%
%=========================================================================================
% Add biogeochemical variables 
if (makenpzd | makepisces | makebioebus | makequota)     makebio = 1;
else     
    makebio = 0;
end
%%
disp(' ')
disp(['Create runoff forcing from Dai and Trenberth''s global monthly climatological run-off dataset'])
disp(' ')
title_name='runoff forcing file (Dai and Trenberth, 2002 dataset)';

%=========================================================================================
define_dir=0 ;  %%->flag to define directly the orientation / direction of the runoff
%
if define_dir==1
    
    %% Define orientation/direction of the flow. First column is the u- (0) or v- (1)
    %% orientation. Second column is the direction of the flow in the choosen orientation
    
% $$$     %%EXAMPLE GIGATl FAMILY for river # 1 to 9 without #7 and #8
% $$$     dir(1,:) = [0 ,  1];  % # Amazon
% $$$     dir(2,:) = [0 , -1];  % # Congo
% $$$     dir(3,:) = [0 ,  1];  % # Orinoco
% $$$     dir(4,:) = [1 , -1];  % # Mississipi
% $$$     dir(5,:) = [0 ,  1];  % # Parana
% $$$     dir(6,:) = [1 ,  1];  % # Tocantins
% $$$     dir(7,:) = [0 ,  1];  % # Tapajos same dir/sense as Amazon
% $$$     dir(8,:) = [0 ,  1];  % # Saint Lawrence
% $$$     dir(9,:) = [0 ,  1];  % # Xingu same dir/sense as Amazon
% $$$     dir(10,:)= [1 ,  1];  % # Magdalena
% $$$     dir(11,:)= [1,  -1];  % # Urugay same dir/sense as Parana
% $$$     dir(12,:)= [0 ,  0];  % # Danube (not used)
% $$$     dir(13,:)= [0 , -1];  % # Niger  

%ARVOR FAMILY for river # 1 to 6 without Danube
    dir(1,:) = [0 ,  1];  % # Orinoco
    dir(2,:) = [0 ,  1];  % # Mississipi
    dir(3,:) = [0 ,  1];  % # Saint Lawrence
    dir(4,:) = [1 ,  1];  % # Magdalena
    dir(5,:) = [0 ,  0];  % # Danube (not used)
    dir(6,:) = [1 , -1];  % # Niger  

end
%=========================================================================================
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[latriv,lonriv,my_flow,myrivername,rivernumber]=runoff_glob_extract(grdname,global_clim_rivername);
%% => rivernumber is the total number of river in the domain

if rivernumber == 0  %%at least a river
    disp(['Create a "fictive" runoff forcing file :  river with no discharge'])
    create_runoff(rivname,grdname,title_name,...
                  qbar_time,qbar_cycle, ...
                  'fictiveriver_@nest',1,18,1,psource_ncfile_ts,makebio,makepisces,makequota)
    my_flow=0;
    nw=netcdf(rivname,'w');
    disp(['Write in runoff file'])
    nw{'Qbar'}(:) = my_flow';
    if psource_ncfile_ts==1
        nw{'temp_src'}(:) = my_temp_src';
        nw{'salt_src'}(:) = my_salt_src';
        if makebio
            nw{'NO3_src'}(:) = my_no3_src';
            if makepisces
               nw{'PO4_src'}(:) = my_po4_src';
               nw{'Si_src'}(:)  = my_sil_src';
               nw{'DIC_src'}(:) = my_dic_src';
               nw{'DOC_src'}(:) = my_doc_src';
               nw{'TALK_src'}(:) = my_talk_src';
               if makequota
                  nw{'DON_src'}(:) = my_don_src';
                  nw{'DOP_src'}(:) = my_dop_src';
               end
            end
        end
    end
    close(nw)
    
    disp([' '])
    disp(['Line to enter in the croco.in file in the psource_ncfile section :'])
    disp(['-----------------------------------------------------------------'])
    disp(['psource_ncfile:   Nsrc  Isrc  Jsrc  Dsrc qbardir  Lsrc  Tsrc   runoff file name'])
    disp(['                           CROCO_FILES/croco_runoff.nc(.#nestlevel)'])
    disp(['                 ',num2str(0)'])
    
    %%
else
    latriv=latriv';
    lonriv=lonriv';
    %%
    rivername=strvcat(myrivername);
    rivernumber=size(rivername,1);
    rivname_StrLen=size(rivername,2);
    %%
    %% Determine the positions of the river, from its lon/lat position
    %% extract j and i to put in croco.in / croco.in.1 for that use of croco_grd.nc
    %%
    %% Read the grid
    %%
    [lat,lon,mask]=read_latlonmask(grdname,'r');
    [latu,lonu,masku]=read_latlonmask(grdname,'u');
    [latv,lonv,maskv]=read_latlonmask(grdname,'v');
    lonmin=min(lon(:));
    lonmax=max(lon(:));
    latmin=min(lat(:));
    latmax=max(lat(:));
    %%
    %% Perform a first pass on each river
    %%
    disp(' ')
    disp(['First guess:'])
    disp(['============'])
    %
    for k= 1:rivernumber
        disp([' '])
        disp(['- Process river #',num2str(k),': ',char(myrivername(k,:))])
        indomain(k)=check_domain_runoff(lon,lat,lonriv(k),latriv(k));
        [j,i]=runoff_grid_pos(lon,lat,lonriv(k),latriv(k));
        J(k)=j-1;
        I(k)=i-1;
        I(I<1)=1;
        J(J<1)=1;

        disp(['  Position is approximetly  J=',num2str(J(k)),' and I=',num2str(I(k))])
        disp(['     lon src in grid (rho point) ~',num2str(lon(J(k),I(k)))])
        disp(['     lat src in grid (rho point) ~',num2str(lat(J(k),I(k)))])
    end
    %==========================================
    mod_manual = 0 ;
    if mod_manual
    % Modify position runoff manually
    % River number to modify number to modify 
        kmod = 3;
        J(kmod) =  J(kmod); 
        I(kmod) = I(kmod) - 10 ; 
        disp(['River #',num2str(kmod),': ',char(myrivername(kmod,:)),' modified manually'])
    end
    %==========================================
    %%
    %% Check the river you really have to process.
    %% Remove rivers out of the domain, if any...
    %%
    rivertoprocess=find(indomain==1);
    %%=================
    %%debug
    %%rivertoprocess=1
    %=================
    number_rivertoprocess=length(rivertoprocess);
    rivername=strvcat(myrivername(rivertoprocess,:));
    rivernumber=number_rivertoprocess;
    rivname_StrLen=size(rivername,2);
 
    %% Make a figure

    if plotting==1
        dl = mean(lat(2:end,1) - lat(1:end-1,1));
        fig1=figure(1);
        set(fig1,'Position',[659 34 1334 1307],'units','pixels');
        m_proj('mercator',...            % dl/2 shift to have cells centered on lon/lat
               'lon',[lonmin-dl/2 lonmax-dl/2],...
               'lat',[latmin-dl/2 latmax-dl/2]);
        m_pcolor(lon-dl/2,lat-dl/2,mask) 
        shading flat
        m_grid('box','fancy','xtick',5,'ytick',5,'tickdir','out');
        set(findobj('tag','m_grid_color'),'facecolor','white');
        hold on
        for k=1:number_rivertoprocess
            lon_src=lon(J(k),I(k));
            lat_src=lat(J(k),I(k));
            [px,py]=m_ll2xy(lon_src,lat_src);
            h1=plot(px,py,'ro');
            set(h1,'Clipping','off')
            legend(h1,'first guess position');
            h2=m_text(lon(J(k),I(k)),lat(J(k),I(k))+0.1,myrivername(k,:));
            set(h2,'fontweight','demi','fontsize',13);
        end
        %pause
        %%
        if plotting_zoom
        %% to be adapted regarding the configuration
            % extendpointY=[500 500 500 500 ...
            %               500 600 900 500 ...
            %               500 500 500 50 ...
            %               500 ];
            % extendpointX=[500 500 500 500 ...
            %               500 600 900 500 ...
            %               500 500 500 50 ...
            %               500 ];
            extendpointY=20*ones(number_rivertoprocess,1);
            extendpointX=20*ones(number_rivertoprocess,1);
            nfig0=5;
            for k=1:number_rivertoprocess
                nfig=nfig0+k;
                figure(nfig)
                
                lon_src=lon(J(k),I(k));
                lat_src=lat(J(k),I(k));
                
                extendpointYkk=extendpointY(k);
                extendpointXkk=extendpointX(k);
                
                Izoom_start=max(1,I(k) - extendpointXkk);
                Izoom_end=min(size(mask,2),I(k) + extendpointXkk); 
                Jzoom_start=max(1,J(k) - extendpointYkk);
                Jzoom_end=min(size(mask,1),J(k) + extendpointYkk); 
                
                %% Generic case    
                mylon_point=lon(Jzoom_start : Jzoom_end,Izoom_start : Izoom_end);
                mylat_point=lat(Jzoom_start : Jzoom_end,Izoom_start : Izoom_end);
                mymask_point=mask(Jzoom_start : Jzoom_end,Izoom_start : Izoom_end);    
                
            %     %% -Replace by specific case for  GIGATL1
            %    if (k ~= 12 )
            %         mylon_point=lon(J(k)+1 - extendpointYkk : J(k)+1 + extendpointYkk,...
            %                          I(k)+1 - extendpointXkk : I(k)+1 + extendpointXkk);
            %         mylat_point=lat(J(k)+1 - extendpointYkk : J(k)+1 + extendpointYkk,...
            %                          I(k)+1 - extendpointXkk : I(k)+1 + extendpointXkk);
            %          mymask_point=mask(J(k)+1 - extendpointYkk : J(k)+1 + extendpointYkk,...
            %                            I(k)+1 - extendpointXkk : I(k)+1 + extendpointXkk);
            %      else
            %          mylon_point=lon(J(k)+1 - extendpointYkk : J(k)+1,...
            %                          I(k)+1 - extendpointXkk : I(k)+1 ) ;
                     
            %          mylat_point=lat(J(k)+1 - extendpointYkk : J(k)+1,...
            %                          I(k)+1 - extendpointXkk : I(k)+1) ;
                     
            %          mymask_point=mask(J(k)+1 - extendpointYkk : J(k)+1,...
            %                               I(k)+1 - extendpointXkk : I(k)+1);
            %      end 
            %      %%- end Specific plotting  
            end
        end
    end
    %%
    %% Choose which river you really want to process...
    %%
    indomain_last=indomain-indomain;
    for k0=1:number_rivertoprocess
        k=rivertoprocess(k0);
        indomain_last(k)=input(['Do you want to use river (Yes[1], No[0]) ?  ', rivername(k0,:)]);
        if indomain_last(k)==0 && plotting_zoom % Close the (zoom) figure if the river is not used
            close(nfig0+k0);
        end
    end
    rivertoprocess=find(indomain_last==1);
    if isempty(rivertoprocess)
        error(['No river selected !'])
    end
    number_rivertoprocess=length(rivertoprocess);
    rivername=strvcat(myrivername(rivertoprocess,:));
    rivname_StrLen=size(rivername,2);
    my_flow=my_flow(:,find(indomain_last==1));
    
    %%
    %% Define the orientation/direction of the flow. First column is the u- (0) or v- (1)
    %% orientation. Second column is the direction of the flow in the choosen orientation
    %%
    
    if define_dir==0
        for k0=1:number_rivertoprocess
            k=rivertoprocess(k0);
            disp(['====='])
            disp(['River ',rivername(k0,:)])
            disp(['Choose the orientation of the flow'])
            dir11=NaN;
            while ~(dir11==0 | dir11==1)
                dir11 = input('0=zonal or 1=meridional. ');
            end
            disp(['Choose the direction of the flow. '])
            
            dir12=NaN;
            while ~(dir12==1 | dir12==-1)
                dir12= input('1 is positive [S-N or W-E], -1 negative [N-S or E-W]. ');
            end
            dir(k,:)=[dir11 dir12];
            disp(['k=',num2str(k)])
            disp(['k0=',num2str(k0)])
            disp(['dir(k,:)=',num2str(dir(k,:))])
        end
    end
    %%
    %% Create the runoff forcing file
    %%
    disp(' ')
    disp(' Create the runoff file...')
    create_runoff(rivname,grdname,title_name,...
                  qbar_time,qbar_cycle, ...
                  rivername,number_rivertoprocess,rivname_StrLen,dir,psource_ncfile_ts,makebio,makepisces,makequota)
    %%
    %% Adjust the rivers positions relative to the mask
    %%
    disp(['Find the real positions of the rivers in the grid: '])
    disp(['==================================================='])
    
    for k0=1:number_rivertoprocess
        k=rivertoprocess(k0);
        disp(['Process final position for river ',myrivername(k,:)])
        disp(['Choose the orientation'])
        jj=J(k);
        ii=I(k);
        dir2=dir(k,:);
        %
        % Compute the u and v cells indexes in CROCO convention
        [jj2for_out,ii2for_out, jj2for,ii2for,jj2,ii2]=locate_runoff(dir2,jj,ii,mask,masku,maskv);
        % 
        J2for_out(k)=jj2for_out;
        I2for_out(k)=ii2for_out;
        %
        J2for(k)=jj2for;
        I2for(k)=ii2for;
        %
        J2(k)=jj2;
        I2(k)=ii2;
        %
        disp([char(myrivername(k,:)),' runoff wet cell (rho-point) indices are :'])
        disp(['matlab             J=',num2str(J2(k)),' and I=',num2str(I2(k))])
        disp(['CROCO fortran conv J=',num2str(J2for(k)),' and I=',num2str(I2for(k))])
        disp([''])
        disp(['For u or v cells indices are (CROCO fortran convention)'])
        disp(['Jsrc=',num2str(J2for_out(k)),' and Isrc=',num2str(I2for_out(k))])
        disp(['===='])
    end
    %%
    %% Adjust the rivers temperature and salinity
    %%
    if psource_ncfile_ts==1
        disp([' '])
        disp(['Adjust the rivers temperature and salinity '])
        %%disp([' Use the closest surface point in the climatology file '])
        my_temp_src0=zeros(rivernumber ,length(woa_time));
        my_salt_src0=zeros(rivernumber ,length(woa_time));
        if makebio==1
            my_no3_src0=zeros(rivernumber,length(woa_time));
            if makepisces==1
               my_po4_src0=zeros(rivernumber,length(woa_time));
               my_sil_src0=zeros(rivernumber,length(woa_time));
               my_dic_src0=zeros(rivernumber,length(woa_time));
               my_doc_src0=zeros(rivernumber,length(woa_time));
               my_talk_src0=zeros(rivernumber,length(woa_time));
               if makequota==1
                  my_don_src0=zeros(rivernumber,length(woa_time));
                  my_dop_src0=zeros(rivernumber,length(woa_time));
               end
            end
        end
        
        if psource_ncfile_ts_auto
            %%==============================================================
            ncclim=netcdf(clmname);
            N=length(ncclim('s_rho'));
            
            for k0=1:number_rivertoprocess
                k=rivertoprocess(k0);
                %%
                %% For temperature, use the target surface rho_point in the clim file
                %% to reduce any heat flux induced by the rivers
                %%
                if dir(k,1)==0 & dir(k,2)==1         %Zonal flux W->E
                   i_rho=I2(k); j_rho=J2(k);
                elseif dir(k,1)==0 & dir(k,2)==-1    %Zonal flux E->W
                   i_rho=I2(k)-1; j_rho=J2(k);                      
                elseif dir(k,1)==1 & dir(k,2)==1     %Meridional flux S->N
                   i_rho=I2(k); j_rho=J2(k);                   
                elseif dir(k,1)==1 & dir(k,2)==-1    %Meridional flux N->S
                   i_rho=I2(k); j_rho=J2(k)-1;                      
                end 
                T=squeeze(ncclim{'temp'}(:,N,j_rho,i_rho));
                my_temp_src0(k,:)=T';
                %%
                %% For salinity ... ?
                %%
                S=squeeze(ncclim{'salt'}(:,N,j_rho,i_rho))-10; % hum...
                S(S<2)=2.0001; % to prevent negative salinities in the
                               % equation of state
                disp(['  Use psource_ncfile_ts_auto using S = sclim -10 '])
                disp(['  Check line 464 in make_runoff.m to change ' ...
                       'this arbitrary runoff salinity'])
                %%S=2;
                my_salt_src0(k,:)=S';
                if makebio==1
                    NO3=squeeze(ncclim{'NO3'}(:,N,j_rho,i_rho));
                    my_no3_src0(k,:)=NO3';
                    if makepisces==1
                       PO4=squeeze(ncclim{'PO4'}(:,N,j_rho,i_rho));
                       my_po4_src0(k,:)=PO4';
                       Si=squeeze(ncclim{'Si'}(:,N,j_rho,i_rho));
                       my_sil_src0(k,:)=Si';
                       DIC=squeeze(ncclim{'DIC'}(:,N,j_rho,i_rho));
                       my_dic_src0(k,:)=DIC';
                       DOC=squeeze(ncclim{'DOC'}(:,N,j_rho,i_rho));
                       my_doc_src0(k,:)=DOC';
                       TALK=squeeze(ncclim{'TALK'}(:,N,j_rho,i_rho));
                       my_talk_src0(k,:)=TALK';
                       if makequota==1
                          DON=squeeze(ncclim{'DON'}(:,N,j_rho,i_rho));
                          my_don_src0(k,:)=DON';
                          DOP=squeeze(ncclim{'DOP'}(:,N,j_rho,i_rho));
                          my_dop_src0(k,:)=DOP';
                       end
                    end
                end
            end
            close(ncclim)
            
            %%==============================================================
        elseif psource_ncfile_ts_manual
            
            %% Alternativaly : Define all mannually the tracer 
            %% t, s, and eventually biogeochemical tracer concentration
            
            %% my_temp/salt/_src0 : line are rivers / rows are month
            %% size(my_temp_src0)=[number_rivertoprocess 12]
            %% my_temp/salt/_src : line are rivers effectivelly processed / rows are month
            
            %% Example for GIGATL FAMILY
            my_temp_src0(:,:)=[28.1677   27.7407   27.7290   27.9662   28.2562 ...
                               28.3312   28.3832   28.4892   28.8444   28.9035 ...
                               28.8618   28.6118 ;
                               
                               27.4397   28.2505   28.7115   28.3359   26.6957 ...
                               24.3065   22.5866   22.1818   23.2263   24.9318 ...
                               26.2492   26.8046 ;
                               
                               26.9818   27.0154   27.2094   27.5274   27.7536 ...
                               27.8428 28.1083   28.5623   28.7645   28.6072 ...
                               28.1333 27.3927 ;
                               
                               18.4467   18.1003   19.3736   21.6955   24.9263   27.6829 ...
                               29.2091   29.5412   28.3085   25.5634   22.6769 ...
                               20.2277 ;
                               
                               22.6310   23.2607   22.4562   20.2833   17.2070   13.8076 ...
                               11.6466   11.2372   12.5740   15.1383   18.1340 ...
                               20.6393 ;
                               
                               27.9732   27.7875   27.9070   28.2647   28.5391   28.5458 ...
                               28.3633   28.1903   28.2218   28.2947   28.3242 ...
                               28.1359 ;
                                   
                               28.1677   27.7407   27.7290   27.9662   28.2562   28.3312 ...
                               28.3832   28.4892   28.8444   28.9035   28.8618 ...
                               28.6118 ;
                               
                               -0.3210   -0.5013   -0.1819    1.5011    5.6784   10.6697 ...
                               14.2692   14.6227   11.4167    7.0775    3.4114 ...
                               1.1195 ;
                               
                               28.1677   27.7407   27.7290   27.9662   28.2562   28.3312 ...
                               28.3832   28.4892   28.8444   28.9035   28.8618 ...
                               28.6118 ;
                               
                               26.9245   26.7369   26.8849   27.2396   27.8931   28.3006 ...
                               28.1773   28.4447   28.7895   28.8083   28.3689 ...
                               27.6033 ;
                               
                               22.6310   23.2607   22.4562   20.2833   17.2070   13.8076 ...
                               11.6466   11.2372   12.5740   15.1383   18.1340 ...
                               20.6393 ;
                               
                               9 9 9 9 9 9 9 9 9 9 9 9 ; %% <= not used
                                   
                               28.6263   28.9582   29.3015   29.4045   29.1542   28.3245 ...
                               27.3215   26.8242   26.8997   27.4027   28.0562 ...
                               28.4404 ];
            
% $$$             my_salt_src0(:,:)=[  31.8610   30.2483   28.0273   25.8570 ...
% $$$                                 25.5184   26.1775   27.5553   28.7191   29.7456 ...
% $$$                                 30.7106   31.5538   32.1978 ;
% $$$                                 
% $$$                                 30.1895   30.2550   30.4903   30.8359   32.1160   33.4259 ...
% $$$                                 33.8053   33.5998   32.8721   32.1840   31.4867 ...
% $$$                                 30.8132 ;
% $$$                                 
% $$$                                 33.1099   33.8338   34.3348   34.4912   34.0427   33.1361 ...
% $$$                                 31.6499   30.3496   30.2760   31.2658   32.1016 ...
% $$$                                 32.5288 ;
% $$$                                 
% $$$                                 31.4802   30.8317   29.5911   28.0779   27.4586   27.6335 ...
% $$$                                 28.4153   29.0401   29.9319   31.0356   31.6389 ...
% $$$                                 31.6811 ;
% $$$                                 
% $$$                                 28.3315   28.3112   28.3833   28.5506   28.9386   29.3127 ...
% $$$                                 29.1937   28.8813   28.5683   28.3198   28.2123 ...
% $$$                                 28.2293 ;
% $$$                                 
% $$$                                 30.6917   28.2332   25.6524   23.7380   24.2345   26.3038 ...
% $$$                                 28.2613   29.4522   30.6267   31.8787   32.6207 ...
% $$$                                 32.4137 ;
% $$$                                 
% $$$                                 31.8610   30.2483   28.0273   25.8570   25.5184   26.1775 ...
% $$$                                 27.5553   28.7191   29.7456   30.7106   31.5538 ...
% $$$                                 32.1978 ;
% $$$                                 
% $$$                                 30.5800   31.1429   31.3964   30.9931   29.8125   28.4991 ...
% $$$                                 27.7825   27.6972   27.8768   28.8042   29.5877 ...
% $$$                                 30.1445 ;
% $$$                                 
% $$$                                 31.8610   30.2483   28.0273   25.8570   25.5184   26.1775 ...
% $$$                                 27.5553   28.7191   29.7456   30.7106   31.5538 ...
% $$$                                 32.1978 ;
% $$$                                 
% $$$                                 35.5176   35.7714   35.9853   35.8987   35.6900   35.5734 ...
% $$$                                 35.7705   35.4552   35.1000   34.6942   34.7754 ...
% $$$                                 35.1856 ;
% $$$                                 
% $$$                                 28.3315   28.3112   28.3833   28.5506   28.9386   29.3127 ...
% $$$                                 29.1937   28.8813   28.5683   28.3198   28.2123 ...
% $$$                                 28.2293 ;
% $$$                              
% $$$                                 33.1350   33.2516   33.2083   33.0608   32.9854   33.0303 ...
% $$$                                 33.3633   33.4677   33.3594   33.2315   33.0206 ...
% $$$                                 33.0312  ];
            %%            
            my_salt_src0(:,:)=ones(rivernumber,12)*5;
            if makebio
% $$$           my_no3_src0=[0 0 0 0 0 0 0 0 0 0 0 0];
            end
        end
        %%
        %% Effectivelly processed rivers
        %%
        my_temp_src= my_temp_src0(find(indomain_last==1),:);
        my_salt_src= my_salt_src0(find(indomain_last==1),:);
        if makebio 
            my_no3_src = my_no3_src0(find(indomain_last==1),:);
            if makepisces
               my_po4_src  = my_po4_src0(find(indomain_last==1),:);
               my_sil_src  = my_sil_src0(find(indomain_last==1),:);
               my_dic_src  = my_dic_src0(find(indomain_last==1),:);
               my_doc_src  = my_doc_src0(find(indomain_last==1),:);
               my_talk_src = my_talk_src0(find(indomain_last==1),:);
               if makequota
                  my_don_src  = my_don_src0(find(indomain_last==1),:);
                  my_dop_src  = my_dop_src0(find(indomain_last==1),:);
               end
            end
        end
    end
    %%
    %%==============================================================
    %% Continue the figure
    
    if plotting==1
        figure(1)
        hold on
        m_proj('mercator',...     %  dl/2 shift to have cells centered on lon/lat
               'lon',[lonmin-dl/2 lonmax-dl/2],...
               'lat',[latmin-dl/2 latmax-dl/2]);
        for k0=1:number_rivertoprocess
            k=rivertoprocess(k0);

            if dir(k,1)==0 %zonal flow (u-grid)
                lon_src=lonu(J2(k),I2for_out(k)); % Adjusted position
                lat_src=latu(J2(k),I2for_out(k)); % + nice plot
            else           %meridional flow (v-grid)
               lon_src=lonv(J2for_out(k),I2(k)); % Adjusted position
               lat_src=latv(J2for_out(k),I2(k)); % + nice plot  
            end
            [px,py]=m_ll2xy(lon_src,lat_src);
            if dir(k,1)==0 & dir(k,2)==1
               h3=plot(px,py,'k>');
            elseif dir(k,1)==0 & dir(k,2)==-1
               h3=plot(px,py,'k<');     
            elseif dir(k,1)==1 & dir(k,2)==1
               h3=plot(px,py,'k^');
            elseif dir(k,1)==1 & dir(k,2)==-1
               h3=plot(px,py,'kv');  
            end    
            set(h3,'Clipping','off')
        end
        legend([h1,h3],{'Approximative first guess river location','final adjusted river location'});
        title({'\bf Location of river in the croco grid';'(from Dai and Trenberth dataset)'});
        
        if plotting_zoom
            nfig0=5;
            for k0=1:number_rivertoprocess
                k=rivertoprocess(k0);
                
                nfig=nfig0+k;
                close(nfig);
                figure(nfig);
                
                lon_src0=lon(J(k),I(k));  % First guess position
                lat_src0=lat(J(k),I(k)); 
                if dir(k,1)==0 %zonal flow (u-grid)
                    lon_src=lonu(J2(k),I2for_out(k)); % Adjusted position
                    lat_src=latu(J2(k),I2for_out(k)); % + nice plot
                else           %meridional flow (v-grid)
                   lon_src=lonv(J2for_out(k),I2(k)); % Adjusted position
                   lat_src=latv(J2for_out(k),I2(k)); % + nice plot  
                end
                %
                lon_src_r=lon(J2(k),I2(k));
                lat_src_r=lat(J2(k),I2(k));
                dlzoom = lat(J2(k),I2(k)) - lat(J2(k)-1,I2(k));
                %
                disp(['Final position of the target wet cell (rho-point) is J=',num2str(J2for_out(k)),' and I=',num2str(I2for_out(k))])
                disp(['   lon src in grid ',num2str(lon_src_r)])
                disp(['   lat src in grid ',num2str(lat_src_r)])  
                disp(['   '])
                %
                extendpointYkk=extendpointY(k);
                extendpointXkk=extendpointX(k);
                %
                Izoom_start=max(1,I(k) - extendpointXkk);
                Izoom_end=min(size(mask,2),I(k) + extendpointXkk); 
                Jzoom_start=max(1,J(k) - extendpointYkk);
                Jzoom_end=min(size(mask,1),J(k) + extendpointYkk); 
                
                %% Generic case for zoom 
                mylon_point=lon(Jzoom_start : Jzoom_end,Izoom_start : Izoom_end);
                mylat_point=lat(Jzoom_start : Jzoom_end,Izoom_start : Izoom_end);
                mymask_point=mask(Jzoom_start : Jzoom_end,Izoom_start : Izoom_end); 
                %%
                m_proj('mercator',...   %  dlzoom/2 shift to have cells centered on lon/lat
                       'lon',[min(min(mylon_point-dlzoom/2)) max(max(mylon_point-dlzoom/2)) ],...
                       'lat',[min(min(mylat_point-dlzoom/2)) max(max(mylat_point-dlzoom/2)) ] );
                
                m_pcolor(mylon_point-dlzoom/2,mylat_point-dlzoom/2,mymask_point) ; %shading flat
                m_grid('box','fancy','xtick',5,'ytick',5,'tickdir','out');
                set(findobj('tag','m_grid_color'),'facecolor','white');
                hold on
                [px0,py0]=m_ll2xy(lon_src0,lat_src0); % First guess position
                h11=plot(px0,py0,'ro');
                [px,py]=m_ll2xy(lon_src,lat_src);     % Adjusted position
                if dir(k,1)==0 & dir(k,2)==1
                    h33=plot(px,py,'k>');
                 elseif dir(k,1)==0 & dir(k,2)==-1
                    h33=plot(px,py,'k<');     
                 elseif dir(k,1)==1 & dir(k,2)==1
                    h33=plot(px,py,'k^');
                 elseif dir(k,1)==1 & dir(k,2)==-1
                    h33=plot(px,py,'kv');  
                 end 
                legend([h11,h33],{'Approximative first guess river location','final adjusted river location'});
                title(['\bf', myrivername(k,:)]);
                %
                % pcolor(mylon_point-dlzoom/2,mylat_point-dlzoom/2,mymask_point) ; %shading flat
                % hold on
                % h11=plot(lon_src0,lat_src0,'ro');
                % if dir(k,1)==0 & dir(k,2)==1
                %     h33=plot(lon_src,lat_src,'k>');
                %  elseif dir(k,1)==0 & dir(k,2)==-1
                %     h33=plot(lon_src,lat_src,'k<');     
                %  elseif dir(k,1)==1 & dir(k,2)==1
                %     h33=plot(lon_src,lat_src,'k^');
                %  elseif dir(k,1)==1 & dir(k,2)==-1
                %     h33=plot(lon_src,lat_src,'kv');  
                %  end 
                % legend([h11,h33],{'Approximative first guess river location','final adjusted river location'});
                % title(['\bf', myrivername(k,:)]);
            end
        end
    end
    
    %% Fill the river discharge and eventually
    %% t/s concentration, no3 concentration
    
    nw=netcdf(rivname,'w');
    disp(['Write in runoff file'])
    %%nw{'runoff_position'}(:,:)=[I2(find(indomain_last==1)) J2(find(indomain_last==1))];
    %%nw{'runoff_direction'}(:,:)=dir(find(indomain_last==1),:);
    %%disp(['... river positions'])
    %%Write qbar, temp,salt and bgc variables conc.
    cff = 1; 
    %
    debugflow=0; 
    if debugflow
        cff=100;
        disp(['================== WARNING ============================'])
        disp(['debugflow is on See line 633 in make_runoff.m          ']) 
        disp(['Use a cff coef to increase the discharge= ',num2str(cff)])
        disp(['======================================================='])
    end
    
    my_flow=cff.*my_flow;
    nw{'Qbar'}(:) = my_flow';

    disp(['... discharges'])
    if psource_ncfile_ts==1
        %% take care : no transpostion needed !!
        nw{'temp_src'}(:) = my_temp_src;
        disp(['... temperature concentration'])
        nw{'salt_src'}(:) = my_salt_src;
        disp(['... salt concentration'])
        if makebio
            nw{'NO3_src'}(:) = my_no3_src';
            disp(['... NO3 concentration'])
            if makepisces
               nw{'PO4_src'}(:) = my_po4_src';
               disp(['... PO4 concentration'])
               nw{'Si_src'}(:) = my_sil_src';
               disp(['... Si concentration'])
               nw{'DIC_src'}(:) = my_dic_src';
               disp(['... DIC concentration'])
               nw{'DOC_src'}(:) = my_doc_src';
               disp(['... DOC concentration'])
               nw{'TALK_src'}(:) = my_talk_src';
               disp(['... TALK concentration'])
               if makequota
                  nw{'DON_src'}(:) = my_don_src';
                  disp(['... DON concentration'])
                  nw{'DOP_src'}(:) = my_dop_src';
                  disp(['... DOP concentration'])
               end 
            end
        end
    end
    if psource_ncfile_ts == 1
        disp([' ==>'])
        disp([' PSOURCE_NCFILE_TS = 1'])
        disp([' Variable river discharge in m3/s +  variable tracer ' ...
              '(temp/salt) concentration '])
        if psource_ncfile_ts_auto
            disp([' PSOURCE_NCFILE_TS_AUTO = 1'])
            disp(['   auto definition of the variable runoff tracer '])
            disp(['   concentration using the nearest point'])
            disp(['   in the climatlogy file'])
        else
            disp([' PSOURCE_NCFILE_TS_MANUAL = 1'])
            disp(['    manual definition of the variable runoff tracer'])
            disp(['    concentration (see example in make_runoff.m directly)'])
        end
        
        disp([' ...Note : '])
        disp([' ... The Tsrc value reported in croco.in are the annual-mean tracer value'])
        disp([' ... It''s just for information !'])
        disp([' ... The Tsrc used are read in the runoff netCDF ' ...
              'file created'])
    else
        disp([' ==> '])
        disp([' PSOURCE_NCFILE_TS = 0'])
        disp([' Variable river discharge in m3/s + constant tracer (temp/salt) concentration ' ...
              'fluxes'])
        disp([' ... The Tsrc value reported in croco.in are the ' ...
              'constant value imposed manually'])
    end
    
    close(nw)
    %%
    %% Line to enter in the croco.in file in the psource section
    %%
    disp([' '])
    disp(['Line to enter in the croco.in file in the psource_ncfile section :'])
    disp(['-----------------------------------------------------------------'])
    disp(['psource_ncfile:   Nsrc  Isrc  Jsrc  Dsrc qbardir  Lsrc  Tsrc   runoff file name'])
    disp(['                           CROCO_FILES/croco_runoff.nc(.#nestlevel)'])
    disp(['                 ',num2str(number_rivertoprocess)])
    for k0=1:number_rivertoprocess
        k=rivertoprocess(k0);
        
        if psource_ncfile_ts==1
            T=mean(my_temp_src0(k,:));
            S=mean(my_salt_src0(k,:));
            if makebio
               N=mean(my_no3_src0(k,:));
               if makepisces
                  P=mean(my_po4_src0(k,:));
                  S=mean(my_sil_src0(k,:));
                  A=mean(my_talk_src0(k,:));
                  Di=mean(my_dic_src0(k,:));
                  Dc=mean(my_doc_src0(k,:));
                  if makequota
                     Dn=mean(my_don_src0(k,:));
                     Dp=mean(my_dop_src0(k,:));
                  end
               end
            end
        else
            T=20;
            S=15;
            if makebio
               N=10;
               if makepisces
                  P=1;
                  S=10;
                  A=2300;
                  Di=2100;
                  Dc=28;
                  if makequota
                     Dn=10;
                     Dp=1;
                  end
               end
            end
        end
        %Print the listing
        if makebio==0
           disp(['                        ',num2str(I2for_out(k)),'  ',num2str(J2for_out(k)),...
                     '  ',num2str(dir(k,1)),'  ',num2str(dir(k,2)),'   T  T   ',...
                     num2str(T),' ',num2str(S)])
        else
           if makepisces==0
              disp(['                     ',num2str(I2for_out(k)),'  ',num2str(J2for_out(k)),...
                     '  ',num2str(dir(k,1)),'  ',num2str(dir(k,2)),'   T  T  T ',...
                     num2str(T),' ',num2str(S),' ',num2str(N)])
           else
              if makequota==0
                 disp(['                  ',num2str(I2for_out(k)),'  ',num2str(J2for_out(k)),...
                        '  ',num2str(dir(k,1)),'  ',num2str(dir(k,2)),' 4*T  2*F  T  F  T  2*F  T  12*F  T  2*F'])
              else
                 disp(['                  ',num2str(I2for_out(k)),'  ',num2str(J2for_out(k)),...
                        '  ',num2str(dir(k,1)),'  ',num2str(dir(k,2)),' 4*T  2*F  T  F  T  2*F  T  12*F  T  2*F  2*T 13*F '])
              end		   
           end
	end
    end
    disp(['-----------------------------------------------------------------'])
    
    %%
    %% Plot the seasonal cycle
    %%
    figure(30)
    if psource_ncfile_ts==1
        subplot(3,2,1)
    end
    hold on
    plot([1:12],my_flow)
    legend(rivername)
    box on, grid on
    title(['\bf Monthly clim of the domain run off'])
    xlabel(['\bf Month']);ylabel(['\bf Discharge in m3/s'])
    set(gca,'Xtick',[0.5:11.5],'XtickLabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D']);
    if psource_ncfile_ts==1
        subplot(3,2,3)
        plot([1:12],my_temp_src)
        box on, grid on
        xlabel(['\bf Month']);ylabel(['\bf Temp [C]'])
        set(gca,'Xtick',[0.5:11.5],'XtickLabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D']);
        %
        subplot(3,2,5)
        plot([1:12],my_salt_src)
        box on, grid on
        xlabel(['\bf Month']);ylabel(['\bf Salt [psu]'])
        set(gca,'Xtick',[0.5:11.5],'XtickLabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D']);
        %
           if makebio
              subplot(3,2,2)
              plot([1:12],my_no3_src)
              box on, grid on
              xlabel(['\bf Month']);ylabel(['\bf NO3 [mmol/m3]'])
              set(gca,'Xtick',[0.5:11.5],'XtickLabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D']);
              %
              if makepisces
                 subplot(3,2,4)
                 plot([1:12],my_doc_src)
                 box on, grid on
                 xlabel(['\bf Month']);ylabel(['\bf DOC [mmol/m3]'])
                 set(gca,'Xtick',[0.5:11.5],'XtickLabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D']);
                 %
                 subplot(3,2,6)
                 plot([1:12],my_sil_src)
                 box on, grid on
                 xlabel(['\bf Month']);ylabel(['\bf Si [mmol/m3]'])
                 set(gca,'Xtick',[0.5:11.5],'XtickLabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D']);
              end
           end
    end

end  %% end at least one river !
