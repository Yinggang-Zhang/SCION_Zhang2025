
time_inso = csvread('Insol-t-0-360ka-day-1-365-lat-(65)-meandaily-La04.txt');
time = time_inso(:,1);
inso = time_inso(:,2);
figure 
hold on
plot(time,inso,'k','displayname',['difference:', num2str(string(max(inso) - min(inso)))]);
[insolation_max,max_id] = findpeaks(inso);
plot(time(max_id),insolation_max,'r','displayname',['difference:', num2str(string(max(insolation_max) - min(insolation_max)))]);
l = legend ;
set(l,'fontsize',8)
set(l,'location','northeast')
xlabel('Time (Ma)')
ylabel('Insolation (w/m2)')







