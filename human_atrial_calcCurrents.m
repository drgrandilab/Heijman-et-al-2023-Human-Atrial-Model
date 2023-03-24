function currents = human_atrial_calcCurrents(t,y,p)

% After running a simulation, feed the time vector and state variables into
% this function to compute ionic currents, etc.
% currents: [I_Na,I_Catot];

currents = [];

for i=1:size(t)
    if ceil(i/1000) == i/1000
        disp(['t = ',num2str(ceil(t(i)))]);
    end
    currents = [currents; human_atrial_model(t(i),y(i,:), p , 'currents')];
end