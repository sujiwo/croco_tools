function  create_runoff(runoffname,grdname,title,...
    qbart,qbarc,rivername,rivernumber,...
    runoffname_StrLen,dir,psource_ncfile_ts,biol,pisces,quota)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Create an empty netcdf runoff file
%       runoffname: name of the runoff file
%       grdname: name of the grid file
%       title: title in the netcdf file
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nw=netcdf(runoffname,'clobber');
%%%result = redef(nw);


%
%  Create dimensions
%

nw('qbar_time') = length(qbart);
nw('n_qbar') = rivernumber;
nw('runoffname_StrLen') = runoffname_StrLen;
nw('one') = 1;
nw('two') = 2;
%
%  Create variables and attributes
%
nw{'qbar_time'} = ncdouble('qbar_time');
nw{'qbar_time'}.long_name = ncchar('runoff time');
nw{'qbar_time'}.units = ncchar('days');
nw{'qbar_time'}.cycle_length = qbarc;

if psource_ncfile_ts
    nw{'temp_src_time'} = ncdouble('qbar_time');
    nw{'temp_src_time'}.long_name = ncchar('runoff time');
    nw{'temp_src_time'}.units = ncchar('days');
    nw{'temp_src_time'}.cycle_length = qbarc;
    
    nw{'salt_src_time'} = ncdouble('qbar_time');
    nw{'salt_src_time'}.long_name = ncchar('runoff time');
    nw{'salt_src_time'}.units = ncchar('days');
    nw{'salt_src_time'}.cycle_length = qbarc;
end
nw{'runoff_name'} = ncchar('n_qbar','runoffname_StrLen');
nw{'runoff_name'}.long_name = ncchar('runoff time');

%% => actually not used/ read in croco.in

%%nw{'runoff_position'} = ncdouble('n_qbar','two');
%%nw{'runoff_position'}.long_name = ncchar('position of the runoff (by line) in the CROCO grid');

%%nw{'runoff_direction'} = ncdouble('n_qbar','two');
%%nw{'runoff_direction'}.long_name = ncchar('direction/sense of the runoff (by line) in the CROCO grid');

nw{'Qbar'} = ncdouble('n_qbar','qbar_time');
nw{'Qbar'}.long_name = ncchar('runoff discharge');
nw{'Qbar'}.units = ncchar('m3.s-1');

if psource_ncfile_ts
    nw{'temp_src'} = ncdouble('n_qbar','qbar_time');
    nw{'temp_src'}.long_name = ncchar('runoff temp conc.');
    nw{'temp_src'}.units = ncchar('deg.celsius');
    
    nw{'salt_src'} = ncdouble('n_qbar','qbar_time');
    nw{'salt_src'}.long_name = ncchar('runoff salt conc.');
    nw{'salt_src'}.units = ncchar('psu');
    
    if biol
        nw{'no3_src_time'} = ncdouble('qbar_time');
        nw{'no3_src_time'}.long_name = ncchar('runoff time');
        nw{'no3_src_time'}.units = ncchar('days');
        nw{'no3_src_time'}.cycle_length = 360;
        
        nw{'NO3_src'} = ncdouble('n_qbar','qbar_time');
        nw{'NO3_src'}.long_name = ncchar('runoff no3 conc.');
        nw{'NO3_src'}.units = ncchar('mmol.m-3');
        if pisces
           nw{'po4_src_time'} = ncdouble('qbar_time');
           nw{'po4_src_time'}.long_name = ncchar('runoff time');
           nw{'po4_src_time'}.units = ncchar('days');
           nw{'po4_src_time'}.cycle_length = 360;

           nw{'PO4_src'} = ncdouble('n_qbar','qbar_time');
           nw{'PO4_src'}.long_name = ncchar('runoff po4 conc.');
           nw{'PO4_src'}.units = ncchar('mmol.m-3');

           nw{'si_src_time'} = ncdouble('qbar_time');
           nw{'si_src_time'}.long_name = ncchar('runoff time');
           nw{'si_src_time'}.units = ncchar('days');
           nw{'si_src_time'}.cycle_length = 360;

           nw{'Si_src'} = ncdouble('n_qbar','qbar_time');
           nw{'Si_src'}.long_name = ncchar('runoff si conc.');
           nw{'Si_src'}.units = ncchar('mmol.m-3');

           nw{'dic_src_time'} = ncdouble('qbar_time');
           nw{'dic_src_time'}.long_name = ncchar('runoff time');
           nw{'dic_src_time'}.units = ncchar('days');
           nw{'dic_src_time'}.cycle_length = 360;

           nw{'DIC_src'} = ncdouble('n_qbar','qbar_time');
           nw{'DIC_src'}.long_name = ncchar('runoff dic conc.');
           nw{'DIC_src'}.units = ncchar('mmol.m-3');

           nw{'doc_src_time'} = ncdouble('qbar_time');
           nw{'doc_src_time'}.long_name = ncchar('runoff time');
           nw{'doc_src_time'}.units = ncchar('days');
           nw{'doc_src_time'}.cycle_length = 360;

           nw{'DOC_src'} = ncdouble('n_qbar','qbar_time');
           nw{'DOC_src'}.long_name = ncchar('runoff doc conc.');
           nw{'DOC_src'}.units = ncchar('mmol.m-3');

           nw{'talk_src_time'} = ncdouble('qbar_time');
           nw{'talk_src_time'}.long_name = ncchar('runoff time');
           nw{'talk_src_time'}.units = ncchar('days');
           nw{'talk_src_time'}.cycle_length = 360;

           nw{'TALK_src'} = ncdouble('n_qbar','qbar_time');
           nw{'TALK_src'}.long_name = ncchar('runoff talk conc.');
           nw{'TALK_src'}.units = ncchar('mmol.m-3');
           if quota
              nw{'don_src_time'} = ncdouble('qbar_time');
              nw{'don_src_time'}.long_name = ncchar('runoff time');
              nw{'don_src_time'}.units = ncchar('days');
              nw{'don_src_time'}.cycle_length = 360;

              nw{'DON_src'} = ncdouble('n_qbar','qbar_time');
              nw{'DON_src'}.long_name = ncchar('runoff don conc.');
              nw{'DON_src'}.units = ncchar('mmol.m-3');
            
              nw{'dop_src_time'} = ncdouble('qbar_time');
              nw{'dop_src_time'}.long_name = ncchar('runoff time');
              nw{'dop_src_time'}.units = ncchar('days');
              nw{'dop_src_time'}.cycle_length = 360;

              nw{'DOP_src'} = ncdouble('n_qbar','qbar_time');
              nw{'DOP_src'}.long_name = ncchar('runoff dop conc.');
              nw{'DOP_src'}.units = ncchar('mmol.m-3');
           end
        end
    end
end
%%%result = endef(nw);

%
% Create global attributes
%
nw.title = ncchar(title);
nw.title = title;
nw.date = ncchar(date);
nw.date = date;
nw.grd_file = ncchar(grdname);
nw.grd_file = grdname;
nw.type = ncchar('CROCO runoff file');
nw.type = 'CROCO runoff file';

%
% Write time variables
nw{'qbar_time'} (:) = qbart;
if psource_ncfile_ts
    nw{'temp_src_time'} (:) = qbart;
    nw{'salt_src_time'} (:) = qbart;
    if biol
         nw{'no3_src_time'} (:) = qbart;
         if pisces
            nw{'po4_src_time'} (:) = qbart;
            nw{'si_src_time'} (:) = qbart;
            nw{'dic_src_time'} (:) = qbart;
            nw{'doc_src_time'} (:) = qbart;
            nw{'talk_src_time'} (:) = qbart;
            if quota
               nw{'don_src_time'} (:) = qbart;
               nw{'dop_src_time'} (:) = qbart;
            end
         end
    end
end
for k=1:rivernumber
    nw{'runoff_name'}(k,:) = rivername(k,:);
end
%
close (nw)
