function y = Ibttle(s,t,k,I)
    if length(s)~=length(t)
        disp('warning: shit ain''t right');
    else
        y = zeros(length(s));
        for j=1:length(s)
            y(j) = integral(@(r) (r.^k).*I(r),s(j),t(j));
        end
    end
        