initialisation;


theta0 = [0;0;0;0];
thetaE = [5;30;0;0.2];

[actU,actW,sols] = linesearchWrapper(uIdx,wIdx,theta0,thetaE,sizes);

figure(1)
hold on
i = 1;
while and(~isempty(sols{i}),i<=length(sols))
    plot(sols{i}.x(1,:),sols{i}.x(2,:),'b')
    i = i + 1;
end



% theta0 = thetaE;
% thetaE = [-12;6;.3;.6];
% 
% uIdx = actU;
% wIdx = actW;
% 
% [actU,actW,sols] = linesearchWrapper(uIdx,wIdx,theta0,thetaE,sizes);
% 
% 
% i = 1;
% while ~isempty(sols{i})
%     plot(sols{i}.x(1,:),sols{i}.x(2,:),'r')
%     i = i + 1;
% end
hold off