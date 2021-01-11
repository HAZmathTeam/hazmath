function stat=plot_lexi(a,x,p,level,lev_plot,num_fig)
    figure(num_fig)
    gplot(a,x); hold on
    axis equal;
    axis off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(lev_plot>0)
        for k1=1:length(p)
            k=p(k1);
            if(level(k)<=lev_plot)
                plot(x(k,1),x(k,2),'ro',x(k,1),x(k,2),'r*');
                title(['level fixed=',int2str(level(k1))]);
            end
        end
    else
        for k1=1:length(p)
            k=p(k1);
            plot(x(k,1),x(k,2),'ro',x(k,1),x(k,2),'r*');
            title(['level=',int2str(level(k1))]);
            pause
        end
    end
    hold off
    stat=0;
end
