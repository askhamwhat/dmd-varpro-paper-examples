
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

set(groot, 'defaultLineMarkerSize',10)
set(groot, 'defaultLineLineWidth',2)

format shorte