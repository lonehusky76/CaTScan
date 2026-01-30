function JM_subplot_inloop_text(time,location,label,figurenumber,p,i,maxcells_per_plot,name)

f = figure(figurenumber);
set(f, 'Name', name);
hold on;

subplot(p(1),p(2),1+mod(i-1,maxcells_per_plot));
text(time,location,label);
hold on;

end

