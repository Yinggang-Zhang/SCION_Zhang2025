%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create time frame and insolation variations for each latitude band %
%%%%%%%%%%%%%%%%%%% by Yinggang Zhang in 02/sep/2021  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load 'latitude_inso_99_89.mat'
load 'INTERPSTACK_2021_improved.mat'

time = -524:0.01:-514;

latitude_inso_99_89 = latitude_inso_99_89';
    
figure
for latitude_id = 1:1:20
    % insolation at each latitude band
    inso_at_lati = latitude_inso_99_89(latitude_id,:);
    % find insolation peaks insolation_max and its corresponding index max_id
    [up1,lo1] = envelope(inso_at_lati,10,'peak');
    %plot insolations over time at each latitude
    subplot(10,2, latitude_id)
    hold on
    box on
    plot(time,inso_at_lati,'k');
    plot(time,up1,'r');
    xlabel('Time (Ma)');
    ylabel([num2str(round(INTERPSTACK.lat(latitude_id),1)), '°'])
end
