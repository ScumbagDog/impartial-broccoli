global zmax;
global numRows;
global numCols;
global SensorReadings;
SensorReadings = Sensor;
[numRows,numCols]=size(SensorReadings);
zmax = 155;
tic;
SensorReadingSize = 500 * 16;
ohit = Ohit();
lambda = LambdaShort();
disp("ohit = "+ohit+" Lambda = "+lambda)
differnceToStop = 0.000000000000001;
oldzhit=0;
oldzshort = 0;
oldzmaxs = 0;
oldzrand = 0;
oldohit = ohit;
oldlambda = lambda;
counter = 0;
while 1 
    ehit = 0;
    eshort = 0;
    emax = 0;
    erand = 0;
    ehitReading=0;
    eshortReading= 0;
    
    for i =1:numCols
        for j= 2:numRows
            PhitCon = Phit(SensorReadings(j,i),SensorReadings(1,i),ohit);
            PshortCon = Pshort(SensorReadings(j,i),SensorReadings(1,i),lambda);
            PrandCon = Prand(SensorReadings(j,i));
            PmaxCon = Pmax(SensorReadings(j,i));
            
            n = (PhitCon+PshortCon+PrandCon+PmaxCon).^(-1);
            
            ehit = ehit + n*PhitCon;
            eshort = eshort + n*PshortCon;
            emax = emax + n*PmaxCon;
            erand = erand + n*PrandCon;
            ehitReading=ehitReading + n*PhitCon*(SensorReadings(j,i)-SensorReadings(1,i).^(-1));
            eshortReading= eshortReading + n*PshortCon*SensorReadings(j,i);
        end
    end 
    disp(" ");
    disp("Round "+counter+":");
    zhit = ehit/SensorReadingSize;
    disp("zhit = "+zhit);
    zshort = eshort/SensorReadingSize;
    disp("zshort = "+zshort);
    zmaxs = emax/SensorReadingSize;
    disp("zmaxs = "+zmaxs);
    zrand = erand/SensorReadingSize;
    disp("zrand = "+zrand);
    disp("IsOne? ="+(zrand+zmaxs+zshort+zhit));
    ohit = sqrt(1/ehit*ehitReading);
    disp("ohit = "+ohit);
    lambda = eshort/eshortReading;
    disp("lambda = "+lambda);
    
    
    shouldStop = ShouldStop(zhit,oldzhit,differnceToStop) && ShouldStop(zshort,oldzshort,differnceToStop) && ShouldStop(zmaxs,oldzmaxs,differnceToStop) && ShouldStop(zrand,oldzrand,differnceToStop)&& ShouldStop(ohit,oldohit,differnceToStop) && ShouldStop(lambda,oldlambda,differnceToStop);
    oldzhit = zhit;
    oldzshort = zshort;
    oldzmaxs = zmaxs;
    oldzrand = zrand;
    oldohit = ohit;
    oldlambda = lambda;
    disp("shouldStop = "+shouldStop);
    if shouldStop
    break;
    end
    counter = counter +1;
    
end

toc;

function o = Ohit()
    global numRows;
    global numCols;
    global SensorReadings;
    ohit_sum = 0;
    ohit_hits = 0;
    for i =1:numCols
        for j= 2:numRows
            if 0 <= SensorReadings(j, i) && SensorReadings(j, i) <= SensorReadings(1, i) 
                ohit_sum = ohit_sum+ (SensorReadings(j, i) - SensorReadings(1, i)).^2;
                ohit_hits = ohit_hits+ 1;
            end
        end
    end 
    o = sqrt(1/ohit_hits*ohit_sum);
end

function l = LambdaShort()
    global numRows;
    global numCols;
    global SensorReadings;
    z_shortSum =0;
    z_hits=0;
    [numRows,numCols]=size(SensorReadings);
    for i=1:numCols
        for j=2:numRows
            z_shortSum =z_shortSum+ SensorReadings(j, i);
            z_hits =z_hits+ 1;
        end
    end
    l=z_hits/z_shortSum;
end

function n = bigN(x,j,z)
    n = exp(-1/2*(x - j).^2/z.^2)/sqrt(2*pi*z.^2);
end
function n = nhit(j,z)
    global zmax;
    n = zmax+j+z;
    %(integral((@(x) exp(-1/2*(x - j).^2/z.^2)/sqrt(2*pi*z.^2)),0,zmax)).^(-1);
end
function p = Phit(x,j,z)
    if 0 <= x && x <= 255
       p = nhit(j,z)*bigN(x,j,z);
    else
        p = 0;
    end
end
function n = nshort(j,lambda)
    n = 1/(1-exp((-lambda)*j));
end
function p = Pshort(x,j,lambda)
    if 0 <= x && x <= j
       p = nshort(j,lambda)*lambda*exp((-lambda)*x); 
    else
    p = 0;
    end
end
function p = Pmax(x)
    global zmax;
    if x >= zmax
        p =1;
    else
    p =0;
    end
end
function p = Prand(x)
    global zmax;
    if 0 <= x && x < zmax
        p = 1/zmax;
    else
    p = 0;
    end
end
function s = ShouldStop(x,y,n)
    if x < y 
        s = (y-x)<=n;
    else 
        s = (x-y)<=n;
    end
end







