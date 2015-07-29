% Data Extraction and plotting function for 2d unsplit nodal DG
% By Devin Light 5/1/14
% ---

function out = plot_2dadv(methname,which_test,ncfilename,res,file_out,stat,contours,axis_lim,print_extrema,saveOutput)
    Qname = strcat('Q',res{1});
    xname = strcat('x',res{1});
    yname = strcat('y',res{1});
    muname = strcat('mu',res{1});
    weightsname = 'qweights';
    nodesname = 'qnodes';
    tname = 'time';
    
    out.data = nc_varget(ncfilename, Qname);
    out.x = nc_varget(ncfilename, xname);
    out.y = nc_varget(ncfilename, yname);
    out.t = nc_varget(ncfilename, tname);
    out.weights = nc_varget(ncfilename,weightsname);
    out.nodes = nc_varget(ncfilename,nodesname);
    out.mu = nc_varget(ncfilename,muname);
    out.N = length(out.nodes) - 1;
    out.method = methname;
    out.test = which_test;
    
    nt = size(out.t,1);
    nx = size(out.x,1);
    
    if(stat==0) % Just output data
    elseif(stat==1) % Make plot at half time and final time
%        scrsz = get(0,'ScreenSize');
%        fig = figure('Position',[1 scrsz(4)/2 3*scrsz(4)/2 scrsz(4)/2]);
        fig = figure();
        
        ax1 = subplot(1,2,1);
        nlvl = 1;%round(nt/2);
        tmp = squeeze(out.data(nlvl,:,:));
        x = out.x; y = out.y;
        contourf(x,y,tmp,contours);
        %colorbar('location','EastOutside')
        axis image;caxis(axis_lim);
        pos = get(ax1,'Position');
        
        if(print_extrema == 1)
            ftitle = ['Max:',num2str(max(tmp(:))),' ; Min:',num2str(min(tmp(:)))];
            title(ftitle,'FontSize',12);
        end
        
        ax2 = subplot(1,2,2);
        tmp = squeeze(out.data(end,:,:));
        contourf(x,y,tmp,contours);
        colorbar('location','EastOutside')
        axis image;caxis(axis_lim);
        
        set(ax2,'Position',[pos(1)+pos(3)+0.05 pos(2) pos(3) pos(4)]);
        pos = get(ax2,'Position');

        B=colorbar;oldpos = get(B,'Position');
        set(B, 'Position', [pos(1)+pos(3)+0.05 0.28 oldpos(3) 1.4*oldpos(4)])
        caxis(axis_lim);

        if(print_extrema == 1)
            ftitle = ['Max:',num2str(max(tmp(:))),' ; Min:',num2str(min(tmp(:)))];
            title(ftitle,'FontSize',12);
        end

        nxny = [num2str(length(out.x)), 'x', num2str(length(out.y))];
        titlemeth = strrep(out.test,'_',' ');
        ftitle = [out.method, ' ; ', titlemeth, ' ; ', 'N=',num2str(out.N), ' ; nx X ny=',nxny];
        suptitle(ftitle);
        
        name = strcat(file_out, res{1},'.pdf');
        if(saveOutput)
            saveas(fig, name, 'pdf');
        end
        
    elseif(stat==2) % Make error plot
        contourf(out.x,out.y,squeeze(out.data(end,:,:))-squeeze(out.data(1,:,:)));
        colorbar('location','EastOutside');
    elseif(stat==3) % Make animation
        scrsz = get(0,'ScreenSize');
        fig = figure('Position',[1 scrsz(4)/2 3*scrsz(4)/4 3*scrsz(4)/4]);
        for i=1:length(out.t)
            tmp = squeeze(out.data(i,:,:)); 
            contourf(out.x,out.y,tmp,contours);
            colorbar('location','EastOutside')
            axis image; caxis(axis_lim);
            ftitle = ['t=',num2str(out.t(i))];title(ftitle);
            pause(0.05);
        end
    end
    

end