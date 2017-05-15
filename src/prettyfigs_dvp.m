
inverted = 0;

bg = 1;
fg = 0;

if (inverted == 1)
    bg = 0;
    fg = 1;
end

set(groot,'defaultAxesColor',[bg bg bg])
set(groot,'defaultFigureColor',[bg bg bg])
set(groot,'defaultAxesYColor',[fg fg fg])
set(groot,'defaultAxesXColor',[fg fg fg])
set(groot,'defaultTextColor',[fg fg fg])
set(groot,'defaultAxesTickLength',[0.04 0.07])

set(groot, 'defaultLineMarkerSize',15)
set(groot, 'defaultLineLineWidth',3)

set(groot, 'defaultScatterLineWidth',1)


format shorte