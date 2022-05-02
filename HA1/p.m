function val = p(x1,x2,y)

    load("stations.mat"); %Maybe not efficient to load each time? Make global
    v=90; gamma=3; sigma=1.5;%Constatns
    %xs=zeros(12,length(x1)); % Placeholders
    xs(1:6,:)=bsxfun(@minus,x1,transpose(pos_vec(1,:)));
    xs(7:12,:)=bsxfun(@minus,x2,transpose(pos_vec(2,:)));
    
%     summ=0;
%     for i=1:6
%         summ=summ+ y(i)-v+10*gamma*log10(vecnorm([x1;x2]-pos_vec(:,i),2,2)).^2;
%     end
%     summ=summ';
    %Y=sum(y);
    %normpdf(y,v-10*gamma*log10(sqrt(xs(1:6,:).^2+xs(7:12,:).^2)),sigma)
    val=prod(normpdf(y,v-10*gamma*log10(sqrt(xs(1:6,:).^2+xs(7:12,:).^2)),sigma));
    %val=1/(2*pi*sigma^2)^3*exp(-1/(2*sigma^2)*summ);
end
