clear all

h = gcf;
axesObjs = get(h, 'Children');
dataObjs = get(axesObjs, 'Children'); 
objTypes = get(dataObjs, 'Type');

xdata = get(dataObjs, 'XData');
ydata = get(dataObjs, 'YData');

y1 = ydata{1};
y2 = ydata{2};

x1 = xdata{1};
x2 = xdata{2};

x1 = x1(12000:end);
x2 = x2(12000:end);
y1 = y1(12000:end);
y2 = y2(12000:end);

f1 = fit(x1',y1',['a*x+b'],'Start',[0 0])
f2 = fit(x2',y2',['a*x+b'],'Start',[0 0])
hold on,plot(x1,y1),plot(x2,y2)

mean([f1.a f2.a])*1e-3

