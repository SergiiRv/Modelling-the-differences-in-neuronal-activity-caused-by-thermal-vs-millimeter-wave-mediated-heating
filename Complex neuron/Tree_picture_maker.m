time = T(2000/0.025+2:end)-2000;
voltage = V(2000/0.025+2:end);
figure('Name','Voltage vs time','NumberTitle','off')
plot(time, voltage, 'k', 'LineWidth', 1);  xlabel('time, msec'); ylabel('Potential, mV');
set(gcf,'color','w'); set(gca,'fontsize',14);%title('Voltage vs time');

%_________________________________________________________________________
A = imread('Voltage vs time');
sizh = size(A,2); % or any other way to find the number of pixels in the horizontal  
                  % direction if you are not using imread.
xmax = sizh*.101;
image([0 xmax],[],A)
%_________________________________________________________________________
hold on
plot([12; 12], [-40; -20], '-k',  [12; 22], [-40; -40], '-k', 'LineWidth', 2)
hold off
axis([[10  34]    -60  140])
text(11.8,-30, '10 mA', 'HorizontalAlignment','right')
text(17,-45, '10 ms', 'HorizontalAlignment','center')
% set(gca, 'XTick', [10:2:34], 'XTickLabel',  {[] 12:2:32 []})    % Used temporarily to get the ‘text’ positions correct
set(gca, 'Visible', 'off')