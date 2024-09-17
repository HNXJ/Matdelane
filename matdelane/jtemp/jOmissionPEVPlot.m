function jOmissionPEVPlot(expvarx, areaX, neuronX, condX, areaNames)

figure;
% areaX = 9;
% neuronX = 53;
% ssX = 2;

a = expvarx{condX(1), areaX};
b = expvarx{condX(2), areaX};
c = expvarx{condX(3), areaX};
d = expvarx{condX(4), areaX};

Mx = max(size(a, 1), size(b, 1));

figure("Position", [0 0 1500 500]);
t = linspace(-970, 5030, 6000);
yt2 = smooth(a(neuronX, :), 200)'*100;
yt2 = yt2 - mean(yt2);
st = ones(1, size(a, 2))/10;
stx = smooth(yt2 + st, 20);
sty = smooth(yt2 - st, 20);
plot(t, stx, "Color", "b", "HandleVisibility", "off");
hold("on");
plot(t, sty, "Color", "b", "HandleVisibility", "off");
patch([t, t(end:-1:1)], [stx;sty(end:-1:1)], [0.5 0.5 1.0], "FaceAlpha", 0.5);
plot(t, yt2, "Color", "k", "HandleVisibility", "off");

yt2 = smooth(b(neuronX, :), 200)'*100;
yt2 = yt2 - mean(yt2);
st = ones(1, size(b, 2))/10;
stx = smooth(yt2 + st, 20);
sty = smooth(yt2 - st, 20);
plot(t, stx, "Color", "r", "HandleVisibility", "off");
hold("on");
plot(t, sty, "Color", "r", "HandleVisibility", "off");
patch([t, t(end:-1:1)], [stx;sty(end:-1:1)], [1.0 0.5 0.5], "FaceAlpha", 0.5);
plot(t, yt2, "Color", "k", "HandleVisibility", "off");

yt2 = smooth(c(neuronX, :), 200)'*100;
yt2 = yt2 - mean(yt2);
st = ones(1, size(c, 2))/10;
stx = smooth(yt2 + st, 20);
sty = smooth(yt2 - st, 20);
plot(t, stx, "Color", "r", "HandleVisibility", "off");
hold("on");
plot(t, sty, "Color", "r", "HandleVisibility", "off");
patch([t, t(end:-1:1)], [stx;sty(end:-1:1)], [0.8 0.7 0.3], "FaceAlpha", 0.5);
plot(t, yt2, "Color", "k", "HandleVisibility", "off");

yt2 = smooth(d(neuronX, :), 200)'*100;
yt2 = yt2 - mean(yt2);
st = ones(1, size(d, 2))/10;
stx = smooth(yt2 + st, 20);
sty = smooth(yt2 - st, 20);
plot(t, stx, "Color", "r", "HandleVisibility", "off");
hold("on");
plot(t, sty, "Color", "r", "HandleVisibility", "off");
patch([t, t(end:-1:1)], [stx;sty(end:-1:1)], [0.6 0.9 0.1], "FaceAlpha", 0.5);
plot(t, yt2, "Color", "k", "HandleVisibility", "off");

title(areaNames{areaX} + " : neuron no. "+ num2str(neuronX));
xlabel("Time (ms)");
ylabel("PEV-change");

xlim([-1000 5000]);
ylim([-0.5 1.5]);

legend("1", "2", "3", "4");

end