function JM_subplot_inloop(time,Ca,x_label,y_label,figurenumber,p,i,markertype,maxcells_per_plot,name)

f = figure(figurenumber);
set(f, 'Name', name);
hold on;

subplot(p(1),p(2),1+mod(i-1,maxcells_per_plot));

if markertype == 1
    plot(time,Ca,'o');
elseif markertype == 2
    plot(time,Ca,'x');
elseif markertype == 3
    plot(time,Ca,'.');
else
    plot(time,Ca);
end
hold on;
xlabel(x_label)
hold on;
ylabel(y_label)
hold on;
title(num2str(i));

end

