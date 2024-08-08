close all; clear all;

filename = 'faktori.txt'; % file with coda Q-factors
delimiter = '\t';
P = importdata(filename,delimiter,1);

for k = 1:2
    j = 0;
    for i = 1:37
        Q1(i) = P.data(1+j, k); % Q-factors for different frequencies
        Q2(i) = P.data(2+j, k);
        Q4(i) = P.data(3+j, k);
        Q6(i) = P.data(4+j, k);
        Q8(i) = P.data(5+j, k);
        Q12(i) = P.data(6+j, k);
        j = j + 6;
    end

    % removing outliers (value > 6*median)
    MED = [median(Q1,'omitnan'), median(Q2,'omitnan'), median(Q4,'omitnan'), median(Q6,'omitnan'), median(Q8,'omitnan'), median(Q12,'omitnan')];

    Q1 = Q1(find(Q1<6*MED(1)));
    Q2 = Q2(find(Q2<6*MED(2)));
    Q4 = Q4(find(Q4<6*MED(3)));
    Q6 = Q6(find(Q6<6*MED(4)));
    Q8 = Q8(find(Q8<6*MED(5)));
    Q12 = Q12(find(Q12<6*MED(6)));

    % removing values outside interquartile range
    LQ = [quantile(Q1,0.25), quantile(Q2,0.25), quantile(Q4,0.25), quantile(Q6,0.25), quantile(Q8,0.25), quantile(Q12,0.25)];
    UQ = [quantile(Q1,0.75), quantile(Q2,0.75), quantile(Q4,0.75), quantile(Q6,0.75), quantile(Q8,0.75), quantile(Q12,0.75)];
    IQR = [iqr(Q1), iqr(Q2), iqr(Q4), iqr(Q6), iqr(Q8), iqr(Q12)];

    Q1 = Q1(find(Q1>=(LQ(1)-1.5*IQR(1))));  Q1 = Q1(find(Q1<=(UQ(1)+1.5*IQR(1))));
    Q2 = Q2(find(Q2>=(LQ(2)-1.5*IQR(2))));  Q2 = Q2(find(Q2<=(UQ(2)+1.5*IQR(2))));
    Q4 = Q4(find(Q4>=(LQ(3)-1.5*IQR(3))));  Q4 = Q4(find(Q4<=(UQ(3)+1.5*IQR(3))));
    Q6 = Q6(find(Q6>=(LQ(4)-1.5*IQR(4))));  Q6 = Q6(find(Q6<=(UQ(4)+1.5*IQR(4))));
    Q8 = Q8(find(Q8>=(LQ(5)-1.5*IQR(5))));  Q8 = Q8(find(Q8<=(UQ(5)+1.5*IQR(5))));
    Q12 = Q12(find(Q12>=(LQ(6)-1.5*IQR(6))));  Q12 = Q12(find(Q12<=(UQ(6)+1.5*IQR(6))));

    % median and mean values
    MED = [median(Q1,'omitnan'), median(Q2,'omitnan'), median(Q4,'omitnan'), median(Q6,'omitnan'), median(Q8,'omitnan'), median(Q12,'omitnan')];
    SR = [mean(Q1,'omitnan'), mean(Q2,'omitnan'), mean(Q4,'omitnan'), mean(Q6,'omitnan'), mean(Q8,'omitnan'), mean(Q12,'omitnan')];

    % median and mean uncertainties
    MSIG = [mad(Q1,1), mad(Q2,1), mad(Q4,1), mad(Q6,1), mad(Q8,1), mad(Q12,1)];
    SSIG = [mad(Q1), mad(Q2), mad(Q4), mad(Q6), mad(Q8), mad(Q12)];
    
    % least square method for median and mean
    fs = [1.5, 3, 6, 9, 12, 18];
    x = log(fs);
    ym = log(MED); ys = log(SR);
    pm = polyfit(x,ym,1); % linear regression for median
    ps = polyfit(x,ys,1); % linear regression for mean
    
    n1m = pm(1); 
    n1s = ps(1); % n value, also b in y = a + b*x
    am = pm(2); 
    as = ps(2); % a in y = a + b*x
    Q01m = exp(am);  
    Q01s = exp(as); % mean Q0

    xks = mean(x)*mean(x); % square of the frequency means
    xsk = mean(x.*x); % mean of frequency squares
    yksm = mean(ym)*mean(ym); 
    ykss = mean(ys)*mean(ys); % square of the median means
    yskm = mean(ym.*ym); 
    ysks = mean(ys.*ys); % mean of median squares

    n1sigm = sqrt(1/6*((yskm-yksm)/(xsk-xks) - n1m^2));
    n1sigs = sqrt(1/6*((ysks-ykss)/(xsk-xks) - n1s^2)); % error for n
    asigm = n1sigm * sqrt(xsk - xks); 
    asigs = n1sigs * sqrt(xsk - xks); % error for a
    Q01sigm = asigm* exp(am); 
    Q01sigs = asigs* exp(as); % error for Q0
    
    % saving values
    nm(k) = n1m; nsigm(k) = n1sigm;
    Q0m(k) = Q01m; Q0sigm(k) = Q01sigm;
    ns(k) = n1s; nsigs(k) = n1sigs;
    Q0s(k) = Q01s; Q0sigs(k) = Q01sigs;
    
end

fileID = fopen('vrijednosti.txt','w'); % writing values to text
fprintf(fileID,'\t%s\t%s\r\n','1. prozor','2. prozor');
fprintf(fileID,'Q0 med\t%.8f\t%.8f\r\n',Q0m);
fprintf(fileID,'Q0 mgr\t%.8f\t%.8f\r\n',Q0sigm);
fprintf(fileID,'n med\t%.8f\t%.8f\r\n',nm);
fprintf(fileID,'n mgr\t%.8f\t%.8f\r\n',nsigm);
fprintf(fileID,'Q0 sr\t%.8f\t%.8f\r\n',Q0s);
fprintf(fileID,'Q0 sgr\t%.8f\t%.8f\r\n',Q0sigs);
fprintf(fileID,'n sr\t%.8f\t%.8f\r\n',ns);
fprintf(fileID,'n sgr\t%.8f\t%.8f\r\n',nsigs);
fclose(fileID);