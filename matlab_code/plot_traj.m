rosinit
exampleHelperROSCreateSampleNetwork

sub = rossubscriber('/robot_traj');
sub2=rossubscriber('/iterations_num');
%pause(1);
msg3=receive(sub2);
j=msg3.Data;

while true

msg3=receive(sub2);
msg2 = receive(sub);

plot(msg2.Data(1),msg2.Data(2),'*','MarkerEdgeColor','r')
hold on
plot(msg2.Data(3),msg2.Data(4),'*','MarkerEdgeColor','g')
hold on
plot(msg2.Data(5),msg2.Data(6),'*','MarkerEdgeColor', 'k')
hold on



end

rosshutdown
