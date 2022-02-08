function drgViolinPlot(x_pos,all_points,edges,rand_offset,which_color_mean,which_color_point,point_size,show_links)
%This does a violin plot for several points

%Do the violin plot
for ii=1:size(all_points,1)
    [mean_out, CIout, x_pos_out(ii,:)]=drgViolinPointPos(all_points(ii,:),edges,x_pos(ii),rand_offset,which_color_mean,which_color_point,point_size);
end

%Show lines linking the points
if show_links==1
    for ii=1:size(all_points,1)-1
        for jj=1:size(all_points,2)
            plot([x_pos_out(ii,jj) x_pos_out(ii+1,jj)],[all_points(ii,jj) all_points(ii+1,jj)],'-','Color',[0.8 0.8 0.8],'LineWidth',2)
        end
    end
end

%Do the violin plot
for ii=1:size(all_points,1)
    [mean_out, CIout]=drgViolinPointPos(all_points(ii,:),edges,x_pos(ii),rand_offset,which_color_mean,which_color_point,point_size,x_pos_out(ii,:));
end