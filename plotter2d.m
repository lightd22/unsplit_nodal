% plotter2d.m
% By: Devin Light
% ------

clc;
clear all;

set(0,'defaultfigureposition',[180 520 2*180 2*180],...
'defaultaxeslinewidth',0.9,'defaultaxesfontsize',8,...
'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.1,...
'defaultlinemarkersize',18); format compact, format long
%%
cd('/Users/Devin/Desktop/R/NodalDG/2d_adv/unsplit_nodal');

tests = {'sine_adv', ... % 1 Uniform adv of sine^4
         'def_cosbell', ... % 2 LeVeque deformation test cosinebell
         'smth_cosbell', ... % 3 Smoother version of LeVeque test
         'fdef_sqwave', ... % 4 LeVeque deformation test square wave
         'hadv_cosinebell', ... % 5 Horizontal advection of steep cosinebell
         'uniform', ... % 6 Uniform field in deformation flow
         'rot_cylinder',... % 7 Solid body rotation of a cylinder
         'rot_cylinder_modified' %8 solid body rotation for comparing to frank's code
         };
res = {'1','2','3','4','5'};

which_test = tests(2);
which_res = res(2);
ncfilename = strcat('dg2d_' ,which_test{1}, '.nc');

contours = -.2:0.1:1.2;axis_lim = [-0.2 1.2];
%contours = 0.0:0.1:2.0;axis_lim = [0.0 2.0];
%% Cycle through methods

stat = 1;
print_extrema = 1;
whichMethods = [3 4 5 6];
subDir = 'n4/';
for nmethod=1:length(whichMethods)
    n = whichMethods(nmethod);
    if(n==1) 
        methname = 'Nodal DG; No Limiting';
        nc = ['_ndgunlim/' subDir ncfilename];
        file = ['figures/nodal/nod_' which_test{1}];
        [nodal_unlim] = plot_2dadv(methname,which_test,nc,which_res,file,stat,contours,axis_lim,print_extrema,0);
    elseif(n==2)
        methname = 'Nodal DG; Zhang and Shu Limiter';
        nc = ['_ndgzhshu/' subDir ncfilename];
        file = ['figures/zshu/zshu_' which_test{1}];
        [nodal_shu] = plot_2dadv(methname,which_test,nc,which_res,file,stat,contours,axis_lim,print_extrema,0);
    elseif(n==3)
        methname = 'Nodal DG; Split; No Limiting';
        nc = ['_ndgsplun/' subDir ncfilename];
        nodal_split = plot_2dadv(methname,which_test,nc,which_res,file,stat,contours,axis_lim,print_extrema,0);
    elseif(n==4)
        methname = 'Nodal DG; Split; Z&S Limiting';
        nc = ['_ndgsplzs/' subDir ncfilename];
        nodal_splZS = plot_2dadv(methname,which_test,nc,which_res,file,stat,contours,axis_lim,print_extrema,0);
    elseif(n==5)
        methname = 'Nodal DG; Split; FCT + TMAR';
        nc = ['_ndgsplfc/' subDir ncfilename];
        nodal_splFC = plot_2dadv(methname,which_test,nc,which_res,file,stat,contours,axis_lim,print_extrema,0);
    elseif(n==6)
        methname = 'Nodal DG; Split; FCT + TMAR';
        nc = ['_ndgspllm/' subDir ncfilename];
        nodal_splLM = plot_2dadv(methname,which_test,nc,which_res,file,stat,contours,axis_lim,print_extrema,0);

    end
end

%% Make convergence plot
clc; close all;
MS = 'MarkerSize';LW = 'LineWidth';FS = 'FontSize';

stat = 0;
print_extrema = 1;

methods = {
    'nodal_split',...
    'nodal_splZS',...
    'nodal_unlim',...
    'nodal_ZS',...
    };
which_test = tests(2);

ncfilename = strcat('dg2d_' ,which_test{1}, '.nc');
nLvls = 1:5;
whichMethods = [2 4];
for istart = 1:length(whichMethods)
    err = [];
    einf = [];
    dxVec = [];
    nmeth = whichMethods(istart);
    methName = methods{nmeth};
    
    for nRun=1:length(nLvls)
        which_res = res(nLvls(nRun));
        if nmeth==1
            methname = 'Nodal; Split; No Limiting';
            nc = ['ndgsplun/', ncfilename];
            file = ['figures/nodal' which_test{1}];   
            out = plot_2dadv(methname,which_test,nc,which_res,file,stat,contours,axis_lim,print_extrema,0);
            meth.(methName).plotStyle = '--b.';
            meth.(methName).MS = 25;
        elseif nmeth==2
            methname = 'Nodal; Split; Z&S Limit';
            nc = ['ndgsplzs/',ncfilename];
            file = ['figures/zshu' which_test{1}];
            out = plot_2dadv(methname,which_test,nc,which_res,file,stat,contours,axis_lim,print_extrema,0);
             meth.(methName).plotStyle = '.-g';
             meth.(methName).MS = 25;
        elseif nmeth==3
            methname = 'Nodal; Unsplit; No Limiting';
            nc = ['ndgunlim/',ncfilename];
            file = ['figures/nodal' which_test{1}];
            out = plot_2dadv(methname,which_test,nc,which_res,file,stat,contours,axis_lim,print_extrema,0);
             meth.(methName).plotStyle = ':rs';
             meth.(methName).MS = 15;
        elseif nmeth==4
            methname = 'Nodal; Unsplit; Z&S Limit';
            nc = ['ndgzhshu/',ncfilename];
            file = ['figures/zshu' which_test{1}];
            out = plot_2dadv(methname,which_test,nc,which_res,file,stat,contours,axis_lim,print_extrema,0);
            meth.(methName).plotStyle = '-cs';
            meth.(methName).MS = 15;
        end
        meth.(methName).label = methname;
        ic = squeeze(out.data(1,:,:));
        final = squeeze(out.data(end,:,:));
            
        %Compute error
        nex = length(ic(:,1))/(out.N+1); dx = 1./nex; 
        ney = length(ic(1,:))/(out.N+1); dy = 1./ney;
        
        qWeights = out.weights; nQuadPts = length(qWeights);
        W = qWeights'*qWeights;
        tmpError = zeros(nex,ney);
        for i=1:nex
            for j=1:ney
                horizPts = 1+(i-1)*nQuadPts:i*nQuadPts;
                vertPts = 1+(j-1)*nQuadPts:j*nQuadPts;
                diff = W.*abs(final(horizPts,vertPts)-ic(horizPts,vertPts)).^2;
                tmpError(i,j) = sum(diff(:));
            end
        end
        nError = sqrt(0.25*dx*dy*sum(tmpError(:)));
        %nError = sqrt(mean( (ic(:)-final(:)).^2 ));
        err = [err nError];
            
        nError = max(abs(ic(:)-final(:)));
        einf = [einf nError];
        dxVec = [dxVec dx];
    end
    
    meth.(methName).error = err;
    meth.(methName).einf = einf;
    meth.(methName).dxVec = dxVec;
end

fig = figure();
hold on,box on;
disp('---');
lstring = {};
for istart = 1:length(whichMethods)
    nmeth = whichMethods(istart);
    methName = methods{nmeth};
    pltNxVals = 1./meth.(methName).dxVec;
    pltErrVals = meth.(methName).error;
    pltStyle = meth.(methName).plotStyle;
    markerSize = meth.(methName).MS;
    
    plot(pltNxVals,pltErrVals,pltStyle,MS,markerSize,LW,1.5);
    lstring = [lstring ;meth.(methName).label];%{['Modal Unlimited'];['Modal PDDG']};
    
    tmp = [methName, ' error vals:', num2str(pltErrVals)];
    disp(tmp);
    appxOrder = polyfit(log10(pltNxVals),log10(pltErrVals),1);
    tmp = ['   Approximate order = ', num2str(appxOrder(1))];
    disp(tmp);
end
disp('---');

% Add reference line
foo = pltNxVals;
m = 4.0;
coeff = (200)*4^(m)*meth.nodal_split.error(1);
refVals = coeff*foo.^(-m);
plot(foo,refVals,'-',LW,1.5,'Color',[0.75 0.75 0.75]);

yVal = (refVals(2)+refVals(3))/2.0;
xVal = foo(2);%(foo(2)+foo(3))/2;
theta = 32;%atan(1.0/m)*(180/(pi));
h = text(xVal,yVal, 'Fourth Order',FS,16);
set(h, 'rotation', -theta);

leg=legend(lstring,'Location','SouthWest');
set(leg,'FontSize',16);

set(gca,FS,12);
set(gca,'XTick',pltNxVals,'XTickLabel',pltNxVals);
xlim([pltNxVals(1)-0.5 pltNxVals(end)+0.5]);
set(gca,'XScale','log','YScale','log');
xlabel('nex',FS,18);
ylabel('L2 Error',FS,18);

%pow = str2double(which_res);
name = strcat('nodalConv2d.pdf');
print(fig,'-dpdf',name);

