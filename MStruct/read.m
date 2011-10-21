%[out,err] = Read (Path,Format,Ratio)
%	Read -	Function reads data
%	Path -	Path and name of dada file
%	Format -	Data format (string) :
%	'dif'		: Difraction data (as for Difpatan)
%				two comment lines, min, max, delta, scale factor, y-data
%   'difc'      : Difraction data (as for Difpatan)
%				two comment lines, min, max, delta, scale factor, y-data.
%               "out" is a cell array that contains in the first cell
%               the first comment line, in the second cell second comment
%               line, in the third 2theta-data and in the fourth
%               cell intensity-data.
%	'xy'		: x,y data in two colloms
%                   the Ratio - parameter has a meaning of a number of
%                   comment lines which will be ignored
%	'difr'	: Difraction data (as for Difpatan)
%				two comment lines, min, max, delta, scale factor, y-data,
%				but only each ratio-th value is considered and added to output
%				data
%   'fit'   : 'Scardi fit' format. Four columns: 1. 2theta 2. measured data
%                                                3. calcul. data 4. difference
%   'raw'   : 'Scardi raw' format.  1. line - comment
%                                   2. line - num. of points  step  origin  end  time
%                                   internsity in one column
%   'dat'   : 'FulProff input' format.  1. line - origin  step  end  comment
%                                       internsity in one column
%   'srf'   : '3D surface plot' format. 1. line - Nx Ny
%                                       x y z
%                                       .....
%   'gnu'   : multuple column data format as for Gnuplot (# - comment)
%   'gnu1'  : same as 'gnu' but comments can be only in the head of the file
%   'matrix': all lines that are not only numbers are ignored
%   'prf'   : FullProf data format ???
%   'brk'   : Bruker FullProf format 
function [out,err] = read (Path,Format,Ratio)
err = '' ;
fid = fopen (Path) ;
if fid==-1
   out=[];
   err = 'Can not open file' ;
   return ;
end
switch Format
case 'difr',
   fgetl(fid) ;
	fgetl(fid) ;
   dd = fscanf (fid, '%f%f%f', [4] ) ;
   y = fscanf (fid, '%f', [1]) ;
   i = 1 ;
   out = y ;
   while (1)
      [y,count] = fscanf (fid, '%f', [1]) ;
      if (count<1)
         break ;
      end ;
      i = i + 1 ;
      if (mod(i,Ratio) == 0)
         out = [out y] ;
      end ;
   end ; 
   out(2,:) = out ;
	out(1,:) = [dd(1):Ratio*dd(3):Ratio*dd(3)*(length(out)-1)+dd(1)] ;
   fclose (fid) ;
case 'dif',
   fgetl(fid) ;
	fgetl(fid) ;
	dd = fscanf (fid, '%f%f%f', [4] ) ;
	out = fscanf (fid, '%f', [1 inf]) ;
	out(2,:) = out/dd(4) ;
	out(1,:) = [dd(1):dd(3):dd(3)*(length(out)-1)+dd(1)] ;
   fclose (fid) ;
case 'difc',
    out{1} = fgetl(fid) ;
	out{2} = fgetl(fid) ;
	dd = fscanf (fid, '%f%f%f', [4] ) ;
	data = fscanf (fid, '%f', [1 inf]) ;
	out{4} = data/dd(4) ;
	out{3} = [dd(1):dd(3):dd(3)*(length(data)-1)+dd(1)] ;
   fclose (fid) ;
case 'xy'
    if exist('Ratio','var')
        for I = 1:Ratio
            fgetl (fid) ;
        end
    end
    out = fscanf (fid, '%f%f\n', [2 inf]) ;
    fclose (fid) ;
case 'fit'
    out = fscanf(fid, '%f%f%f%f', [4 inf]) ;
    fclose (fid) ;
case 'raw',
   fgetl(fid) ;
	dd = fscanf (fid, '%d%f%f%f%f', [5] ) ;
	out = fscanf (fid, '%f', [1 inf]) ;
	out(2,:) = out/dd(5) ;
	out(1,:) = [dd(3):dd(2):dd(2)*(dd(1)-1)+dd(3)] ;
   fclose (fid) ;
case 'dat',
   line = fgetl(fid) ;
    dd = sscanf (line, '%f%f%f', [3] ) ;
	out = fscanf (fid, '%f', [1 inf]) ;
	out(2,:) = out ;
	out(1,:) = [dd(1):dd(2):dd(3)] ;
   fclose (fid) ;
case 'srf'
	nx = fscanf (fid, '%d', [1] ) ;
    ny = fscanf (fid, '%d', [1] ) ;
	data = fscanf (fid, '%f%f%f', [3 inf]) ;
	out{1} = reshape(data(1,:)',ny,nx);
    out{2} = reshape(data(2,:)',ny,nx);
    out{3} = reshape(data(3,:)',ny,nx);
   fclose (fid) ;
case 'gnu',
    out = [];
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end    
        if (length(tline)<1 | tline(1)=='#'), continue, end,
        out(end+1,:) = str2num(tline);
    end
    out = out.';
    fclose(fid);
case 'gnu1',
    out = [];
    % read comments and empty lines in the file head
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end    
        if isempty(tline), continue, end
        if tline(1)=='#', continue, else break, end
    end
    % read data from the first data line
    if ischar(tline) && ~isempty(tline), out(:,1) = str2num(tline).'; end
    % read all remaining data
    t = fscanf(fid,'%f',[size(out,1),inf]);
    out = [out, t];
    fclose(fid);
case 'matrix',
    out = [];
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end,
        d = str2num(tline);
        if isempty(d), continue, end,
        if (size(d,1)>1), d = d.'; end,
        out(end+1,:) = d;
    end
    out = out.';
    fclose(fid);
case 'prf',
    out = [];
    % comment line
    tline =  fgetl(fid);
    % ??? ??? zero_shift
    tline =  fgetl(fid);
    % 2Tend 2Tstart 2Tstep ??? ???
    tline =  fgetl(fid);
    d = sscanf(tline,'%f',[1 3]);
    t2_end = d(1); t2_start = d(2); t2_step = d(3);
    n = round((t2_end-t2_start)/t2_step)+1;
    % nb.of.phases nb.of.points Lambda1 Lambda2 zero_shift ??? ???
    tline =  fgetl(fid); d = str2num(tline);
    if (d(2)~=n), warning('incompatible nb. of points'), end,
    % nb.diffractions_phase1 nb.diffractions_phase2 ... ??? ???
    line =  fgetl(fid);
    % Iobs
    d = []; d(1,:) = fscanf(fid,'%f',[1 n]);
    % Icalc
    d(2,:) = fscanf(fid,'%f',[1 n]);
    % large integer number - same count as all diffractions
    % positions of peaks (real n.) - same count as all diffractions
    %
    out = [t2_start+[0:n-1]*t2_step; d];
    %
    fclose(fid);
case 'brk',
    out = {};
    dset = [];
    xinfo = [];
    % read lines until the end of file
    nline = 0;
    nlinexinfo = 0;
    while 1,
        tline = fgetl(fid);
        nline = nline + 1;
    	if ~ischar(tline), break, end
        tline = strtrim(tline);
        tline = deblank(tline);
        if isempty(tline), continue, end
        % analyse line
        d = sscanf(tline,'%f',[1 inf]);
        if numel(d)>1,
            % new data set
            % save the last data set
            if numel(dset)>=1,
                x = (0:length(dset)-1)*xinfo(2) + xinfo(1);
                if abs(x(end)-xinfo(3))>1e-6,
                    str = sprintf('line: %d, scan num. %d: incompatible x-values.',nlinexinfo,length(out)+1);
                    warning(str),
                end
                out{end+1} = [x; dset].';
            end
            % get the new header
            if numel(d)~=3,
                error('line: %d, scan num. %d: wrong format of x-info',nline,length(out)+1);
            else
                xinfo = d;
                nlinexinfo = nline;
                dset = [];
            end
        else
            % mayby data value
            if numel(d)~=1,
               error('line %d, scan num. %d: unexpected data format.',nline,length(out)+1);
            else
                % really data value
                dset(end+1) = d;
            end
        end
    end
    % save the last data set
    if numel(dset)>=1,
        x = (0:length(dset)-1)*xinfo(2) + xinfo(1);
        if abs(x(end)-xinfo(3))>1e-6,
            str = sprintf('line: %d, scan num. %d: incompatible x-values.',nlinexinfo,length(out)+1);
            warning(str),
        end
        out{end+1} = [x; dset].';
    end
    % if only one data set - make output more compact
    if numel(out)==0,
        out = [];
    elseif numel(out)==1,
        out = out{1};
    end
otherwise,
   err = 'Bad Format parameter' ; 
end
out = out.' ;


