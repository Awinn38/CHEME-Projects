%% n-anova test
CoffeeA = xlsread('CoffeeLab.xlsx'); % This reads collected data from excel file
CoffeeA(:,4) = []; % Removes the unwanted 4th column 
CoffeeA_R = reshape(CoffeeA,[1,24]); % This reshapes the created matrix into a row vector 
Concentration = { '1Cup' '2Cup' '1Cup' '2Cup' '1Cup' '2Cup' '1Cup' '2Cup'}; % The denotation of levels between concentration
Gradulesize = { 'Course' 'Course' 'Course' 'Course' 'Fine' 'Fine' 'Fine' 'Fine'}; % The denotation of levels between gradule size
Bean_Type = { 'Original' 'Original' 'Dark' 'Dark' 'Original' 'Original' 'Dark' 'Dark'}; % The denotation of levels between bean type
Concentration = [ Concentration Concentration Concentration ]; % Denotation of amount of trials
Gradulesize  = [ Gradulesize Gradulesize Gradulesize]; % Denotation of amount of trials
Bean_Type = [ Bean_Type Bean_Type Bean_Type ]; % Denotation of amount of trails
[~,~,stats] = anovan(CoffeeA_R,{Concentration,Gradulesize,Bean_Type},'model','full','varnames',{'Concentration','Gradule_Size','Bean_Type'}); %n-way anova
subplot(2,1,1)
multcompare(stats,'Dimension',[1 3],'Ctype','bonferroni'); % Bonferroni multiple comparison
%% Absorbance Scatter Plot Original
TrialsA = linspace(1,8,8); % Setting x-axis for general trail count
TrialsO = linspace(1,4,4); % Setting x-axis for Original absorbance trials
AbsorbanceO = zeros(1,4)
for jj = 1:4
        if jj <= 2 
            AbsorbanceO(jj) = mean(CoffeeA(jj,:))
        else if jj > 2
                AbsorbanceO(jj) = mean(CoffeeA(jj+4,:))
            end
        end
end   
p = polyfit(TrialsO,AbsorbanceO,1);
AbsorbancePolyO = polyval(p,TrialsO);
scatter(TrialsO,AbsorbanceO) % Scatter plot of absorbance
title('Absorbance of CGA Within Original Brand')
ylabel('Mean Absorbance(nm)')
hold on 
plot(TrialsO,AbsorbancePolyO) % Linear regression of data
str = {'Condition 1' 'Condtion 2' 'Condition 5' 'Condition 6'};
text(TrialsO,AbsorbanceO,str); % Naming of Data points
hold off
%% Absorbance Scatter plot Dark
TrialsD = linspace(1,4,4); % Setting x-axis for Original absorbance trials
AbsorbanceD = zeros(1,4)
for jj = 1:4
        if jj <= 2 
            AbsorbanceD(jj) = mean(CoffeeA(jj+2,:))
        else if jj > 2
                AbsorbanceD(jj) = mean(CoffeeA(jj+4,:))
            end
        end
end   
p = polyfit(TrialsD,AbsorbanceD,1);
AbsorbancePolyD = polyval(p,TrialsD);
scatter(TrialsD,AbsorbanceD) % Scatter plot of absorbance
title('Absorbance of CGA Within Dark Brand')
ylabel('Mean Absorbance(nm)')
hold on 
plot(TrialsD,AbsorbancePolyD) % Linear regression of data
str = {'Condition 3' 'Condtion 4' 'Condition 7' 'Condition 8'};
text(TrialsD,AbsorbanceD,str); % Naming of Data points
hold off
%% Absorbance Bar Graph
Absorbance = zeros(1,8); % Setting of zeros vector
for ii= 1:8
TrialsA = linspace(1,8,8); % Setting x-axis for trials 
Absorbance(ii)= mean(CoffeeA_R(:,ii)); % Filling of zeros vector
end
AbsorbanceErr = zeros(1,8)
for kk = 1:8
    AbsorbanceErr(kk) = std(CoffeeA(kk,:));
end
bar(Absorbance,'w')
%title('Absorbance of CGA Within Sample')
xlabel('Conditions')
ylabel('Mean Absorbance(nm)')
hold on
%e=std(AbsorbanceErr)*ones(size(TrialsA))
e=std(CoffeeA')
errorbar(Absorbance,e,'o','Color','r')
%% Concentration anova
CoffeeC = CoffeeA_R./(27025*.10)
Concentration = { '1Cup' '2Cup' '1Cup' '2Cup' '1Cup' '2Cup' '1Cup' '2Cup'}; % The denotation of levels between concentration
Gradulesize = { 'Course' 'Course' 'Course' 'Course' 'Fine' 'Fine' 'Fine' 'Fine'}; % The denotation of levels between gradule size
Bean_Type = { 'Original' 'Original' 'Dark' 'Dark' 'Original' 'Original' 'Dark' 'Dark'}; % The denotation of levels between bean type
Concentration = [ Concentration Concentration Concentration ]; % Denotation of amount of trials
Gradulesize  = [ Gradulesize Gradulesize Gradulesize]; % Denotation of amount of trials
Bean_Type = [ Bean_Type Bean_Type Bean_Type ]; % Denotation of amount of trails
[~,~,stats] = anovan(CoffeeC,{Concentration,Gradulesize,Bean_Type},'model','full','varnames',{'Concentration','Gradule_Size','Bean_Type'}); %n-way anova
subplot(2,1,1)
multcompare(stats,'Dimension',[1 3],'Ctype','bonferroni'); % Bonferroni multiple comparison
%% Standard Deviation Plot 
scatter(TrialsA,AbsorbanceErr)
xlabel('Conditions')
ylabel('Absorbance Deviation')