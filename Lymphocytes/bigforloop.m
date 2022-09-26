%% Code for Big for loop
tic
"Starting Big Loop"
% use rats(3), 1, rats(11), rats(13)
% figure(11);clf
% subplot(2,2,4)
% SF2_e = SF2_t / 1;
% parfor m = 1:size(I0_orig,1)
% parfor m = 1:25
% for m = 38
    parfor m = 1
        "Point " + m
        start_time = 0;
%         I0 = I0_orig(m, :); % initial conditions
%         SFends{m, 1} = I0_orig(m, :)
%         ronr_a(m, 1:2) = I0_orig(m, :);
%         ronr_a(m, 4) = LocalFailure(m);
        dose = dose_vec(m);
        num_fx = num_fx_vec(m);
        % use 1, rats(4), rats(11), rats(13)
%         SF2_e = SF2_t(m) / 1;
%         prtk_t = 1-SF2_t(m);
%         prtk_e = 1-SF2_e(m);
%         subplot(1,1,1)
        % plotting the initial conditions
%         x_I0 = I0(1)
%         y_I0 = I0(2)
%         x_I0 = I0(1)/max(x)*(Npoints-1);
%         y_I0 = I0(2)/max(y)*(Npoints-1);
%         plot(x_I0, y_I0,'b.','linewidth',1.5,'markersize',10)
%         xlim([1e-4, 1e8])
%         ylim([1e-10, 1.5e5])
%         hold on
%         if LocalFailure(m) == 0
%             plot(x_I0, y_I0, 'bo', markersize = bubbleSize)
%         end
    
%         hold on
        xpost = 0;
        ypost = 0;
        xpath = 0;
        ypath = 0;
        rnr = 0;
%         SF2_e = SF2_t/rstar(m)
        for rt = 1:size(SFrats, 2)
            q = SFrats(rt);
            SF2_e = SF2_t/q;
            I0 = I0_orig(m, :); % initial conditions
            x_I0 = I0(1)/max(x)*(Npoints-1);
            y_I0 = I0(2)/max(y)*(Npoints-1);
        % how many fractions to run
            for k = 1:num_fx
%                 "Num_fx = " + k
                % recalculate initCond
                initCond = [I0(1)*SF2_e(m) I0(2)*SF2_t(m)];
                sols = solve(initCond);
                y_vec = sols.y(2,:)/max(y)*(Npoints-1);
                x_vec = sols.y(1,:)/max(x)*(Npoints-1);
               if k < num_fx % is this the final fraction?
                   % updating I0 to treat w/ another dose 
                   I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
                   start_time = start_time + fx_dt;
    %                plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'b+','linewidth',0.5,'markersize',2)
    %                hold on
    %                plot(x_vec(1:fx_dt),y_vec(1:fx_dt),'k','linewidth',1); 
               elseif k == num_fx % what to do on final fraction
                   % run it out all the way
    %                plot(x_vec,y_vec,'r-','linewidth',0.1);
                   x_rt = initCond(1)/max(x)*(Npoints-1);
                   y_rt = initCond(2)/max(y)*(Npoints-1);
                   xpost = x_rt;
                   ypost = y_rt;
                   xpath = x_vec; 
                   ypath = y_vec;
    %                plot(x_rt, y_rt ,'r.','linewidth',20,'markersize',10)
    %                if LocalFailure(m) == 0
    %                    plot(x_rt, y_rt, 'ro', markersize = bubbleSize)
    %                end
    %                plot([x_rt, x_I0], [y_rt, y_I0],  "b-", LineWidth=0.1, LineStyle="--")
    %                hold on
                   start_time = start_time + fx_dt;
                   
                   % Percent classification
%                    if y_vec(end) > 10 % y value for resolution
%                        ronr_a(m, 3) = 1;
%                    end
               end
            end
            postRT(m, :) = [xpost, ypost];
            xpaths{m} = xpath;
            ypaths{m} = ypath;
            if ypath
                if ypath(end) > threshold
                     rnr = 1; % represents LRF
                end
            end
            SFends{m, 3}(rt) = rnr
        end
    end
% set(gca, 'xscale', 'log')
% set(gca, 'yscale', 'log')
% xlabel("E", FontSize = 15)
% ylabel("T", FontSize = 15)
% xlim([1e-4, 1e8])
% ylim([1e-10, 1.5e5])
% % tit = "Effects of radiation therapy, rho = " + c3_rho_star
% tit = "SF2_E = SF2_T / " + rats(13);
% title(tit, FontSize = 20)
toc
