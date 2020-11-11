%***************************************************************************************************%

clear all;
close all;
clc;

%***************************************************************************************************%

green = [0 204 0]/255;
grey = [0 0 0]+0.25;

%***************************************************************************************************%

fichier = load('-ascii', 'gain.txt');

freq = fichier(:,1);
gain = fichier(:,2);

lbl_dwn = .1*max(gain);

%************************************ Definition des gabarits *************************************%


%Rectange de gauche
x1_rect1 = 7;
y1_rect1 = -15;
x2_rect1 = 100;
y2_rect1 = 20;
w_rect1 = x2_rect1-x1_rect1;
h_rect1 = y2_rect1-y1_rect1;

%Rectange du milieu
x1_rect2 = 2000;
y1_rect2 = -60;
x2_rect2 = 7000;
y2_rect2 = 5;
w_rect2 = x2_rect2-x1_rect2;
h_rect2 = y2_rect2-y1_rect2;

%Rectange de droite
x1_rect3 = 2*10^5;
y1_rect3 = -20;
x2_rect3 = 2*10^6;
y2_rect3 = 20;
w_rect3 = x2_rect3-x1_rect3;
h_rect3 = y2_rect3-y1_rect3;


%**************************************** Plot en ligne ********************************************%


figure(1)

semilogx (freq, gain, 'color',green, 'linewidth', 1)
axis ([ x1_rect1 x2_rect3 min(gain)-5 max(gain)+5 ])

%axis auto
grid on
hold on

%letitsnow(figure(1))

title ('La reponse en frequences d''un filtre passe-bande du deuxieme ordre','FontSize',15,'FontWeight','bold','Color','b')
xlabel('f (Hz)','FontSize',15,'FontWeight','bold','Color','b')
ylabel ('G(dB)','FontSize',15,'FontWeight','bold','Color','b')
legend ('G_d_B (f)')


%****************************************** Plot en points *****************************************%


%for i=1:length(freq)
%  if freq(i)<100
%    if gain(i)<-15
%      semilogx (freq(i), gain(i), 'o', 'markersize', 6, 'MarkerEdgeColor', green, 'MarkerFaceColor', green)
%    else semilogx (freq(i), gain(i), 'r*', 'markersize', 6)
%    end
%  elseif freq(i)>2000 && freq(i)<7000
%    if gain(i)>5
%      semilogx (freq(i), gain(i), 'o', 'markersize', 6, 'MarkerEdgeColor', green, 'MarkerFaceColor', green)
%    else semilogx (freq(i), gain(i), 'r*', 'markersize', 6)
%    end
%  elseif freq(i)>200000
%    if gain(i)<-20
%      semilogx (freq(i), gain(i), 'o', 'markersize', 6, 'MarkerEdgeColor', green, 'MarkerFaceColor', green)
%    else semilogx (freq(i), gain(i), 'r*', 'markersize', 6)
%    end
%  else
%    semilogx (freq(i), gain(i), 'o', 'markersize', 6, 'MarkerEdgeColor', green, 'MarkerFaceColor', green)
%  end
%  if freq(i)==100 || freq(i)==2000 || freq(i)==7000 || freq(i)==200000
%% (freq(i)<110&&freq(i)>90) || (freq(i)<2500&&freq(i)>1800) || (freq(i)<7500&&freq(i)>6800) || (freq(i)<250000&&freq(i)>180000)  
%    label = num2str(gain(i)); 
%    text (freq(i), gain(i), label, 'Position', [freq(i) gain(i)]);
%  end
%  hold on
%  grid on
%end


%********************************* Affichage gabarit en gris ***************************************%


%clf;
%hold on

rect1 = rectangle ("Position", [x1_rect1, y1_rect1, w_rect1, h_rect1], "FaceColor", grey, "EdgeColor", 'none');
rect2 = rectangle ("Position", [x1_rect2, y1_rect2, w_rect2, h_rect2], "FaceColor", grey, "EdgeColor", 'none');
rect3 = rectangle ("Position", [x1_rect3, y1_rect3, w_rect3, h_rect3], "FaceColor", grey, "EdgeColor", 'none');

set( get(rect1, 'children'), 'facealpha', 0.25 );
set( get(rect2, 'children'), 'facealpha', 0.25 );
set( get(rect3, 'children'), 'facealpha', 0.25 );


%***************************************************************************************************%


saveas (figure(1), 'bode_graph', 'png');
%print -dsvg figure1.svg
%print -dpngalpha myplot.png




