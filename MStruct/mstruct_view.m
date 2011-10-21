function [] = mstruct_view(filename,Ihkl_file_list,xl)

% check input parameters
if ~exist('xl','var'), xl = []; end

if ~exist('Ihkl_file_list','var'), Ihkl_file_list = {}; end

if ~exist('bplotindexes','var') || isempty(bplotindexes) ...
        && ~empty(Ihkl_file_list), bplotindexes=1; end,

bplotdiff = 1;
%filename = 'TiO2-1.54_350-15';

if ~exist('bplotdiff','var') || isempty(bplotdiff), bplotdiff=0; end,
d=read([filename '.dat'],'gnu1');

hold off,

x = d(:,1).';
y = d(:,2).';

% plot data limits
if isempty(xl), xl = [x(1) x(end)]; end
Lindx = x>= xl(1) & x<=xl(2);

% plot data
plot(d(Lindx,1),sqrt(d(Lindx,2)),'ro'),

hold on,

plot(d(Lindx,1),sqrt(d(Lindx,3)),'b','LineWidth',2),

if logical(bplotdiff),
    plot(d(Lindx,1),sign(d(Lindx,4)).*sqrt(abs(d(Lindx,4))),'g','LineWidth',1),
end,

%xlim([d(1,1) d(end,1)])
%xl = xlim;

% background
if (0)
db=read([filename '.bgr'],'xy'); 
%db = read('el_bkg.bgr','xy');
plot(db(:,1),sqrt(db(:,2)),'c'),plot(db(:,1),sqrt(db(:,2)),'co')
end

if(0)
xlim(xl);

set(gca,'FontSize',12)

YTick = get(gca,'YTick');
YTick1 = sign(YTick).*YTick.^2;
YTickLabel = [];
for m=1:length(YTick1),
    YTickLabel = strvcat(YTickLabel,num2str(YTick1(m)));
end,
set(gca,'YTickLabel',YTickLabel),
end

% difraction indexes

cols = {'k','m','c','y','g','r','b'};

% plot diffraction lines indexes
if bplotindexes,
    
    % for each item in Ihkl_file_list print diffractions indexes
    for ifrac=1:numel(Ihkl_file_list),

        fid = fopen(['Ihkl_' Ihkl_file_list{ifrac} '.dat'],'rt'); n = 0; str = '';
        d1 = []; n = 0;
        while logical(1),
            tline = fgetl(fid);
            if isempty(tline) | (tline(1)=='#' & n>0), break, end,
            if (tline(1)=='#'), continue, end,
            n = n + 1; d1(n,:) = str2num(tline);
        end,
        fclose(fid);
        
        Lind = d1(:,5)>1e-7 & d1(:,4)>=xl(1) & d1(:,4)<=xl(2);
        hkl = d1(Lind,1:3); th2 = d1(Lind,4); % fhkl2 = d1(Lind,5)
        yhkl = 1.2*interp1(d(:,1),d(:,2),th2);
        x=[]; y = []; shkl = {};
        for n=1:length(th2),
            x(n) = th2(n); y(n) = yhkl(n);
            shkl{n} = sprintf('- %d %d %d',hkl(n,:));
        end,

        text(x,sqrt(y),shkl,'Rotation',90,'Color',cols{mod(ifrac,numel(cols)+1)}),

    end % ifrac

end % bplotindexes

% axes format and labels
set(gca,'FontSize',12)

xlabel('2\Theta (deg)','FontName','Times','FontSize',16),
ylabel('Intensity (counts)','FontName','Times','FontSize',16),

% axes limits
dx = mean(diff(d(:,1)));
xl = [xl(1)-10*dx xl(end)+10*dx];
xlim(xl)

% set correct axis tick labels
yl = ylim;
set(gca,'FontName','Times','FontSize',12)
ylim(yl);

YTick = get(gca,'YTick');
set(gca,'YTickLabelMode','manual')
set(gca,'YTickMode','manual')
YTick1 = sign(YTick).*YTick.^2;
YTickLabel = [];
for m=1:length(YTick1),
    YTickLabel = strvcat(YTickLabel,num2str(YTick1(m))); %#ok<VCAT>
end,
set(gca,'YTickLabel',YTickLabel),

% title
title(filename,'Interpreter','none','FontName','Times','FontSize',16);

return,

d=read('profileP.dat','gnu');
plot(d(:,1),sqrt(400+250*d(:,2)),'m*')

return,