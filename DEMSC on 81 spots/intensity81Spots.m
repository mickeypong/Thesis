% 1. Read and show the fluorescent image of 81 dots
fileAbrs{1} = 'A2 Amine G25-81strands Gfilter 2s';
fileAbrs{2} = 'D2 Amine NG25-F1-NF234 again Gfilter 2s';
fileAbrs{3} = 'D2 Amine NG25-F1-NF234 again Ofilter 2s';
fileAbrs{4} = 'D2 Amine NG25-F1-NF234 again Rfilter 2s';
fileAbrs{5} = 'D2 Amine NG25-F1-NF234 again Yfilter 1s';
fileAbrs{6} = 'E2 Amine F2-NF134 Gfilter 2s';
fileAbrs{7} = 'E2 Amine F2-NF134 Ofilter 2s';
fileAbrs{8} = 'E2 Amine F2-NF134 Rfilter 2s';
fileAbrs{9} = 'E2 Amine F2-NF134 Yfilter 1s';
fileAbrs{10} = 'F2 Amine F3-NF124 Gfilter 2s';
fileAbrs{11} = 'F2 Amine F3-NF124 Ofilter 2s';
fileAbrs{12} = 'F2 Amine F3-NF124 Rfilter 2s';
fileAbrs{13} = 'F2 Amine F3-NF124 Yfilter 1s';
fileAbrs{14} = 'G2 Amine F4-NF123 Gfilter 2s';
fileAbrs{15} = 'G2 Amine F4-NF123 Ofilter 2s';
fileAbrs{16} = 'G2 Amine F4-NF123 Rfilter 2s';
fileAbrs{17} = 'G2 Amine F4-NF123 Yfilter 1s';

for j = 1:17
    fileAbr = fileAbrs{j};
    imgFilename = strcat('/Users/jackie/Downloads/SGTC microscope/' ...
        ,'2016-07-22-81 strands Atto sequence color/' ...
        ,'registered grey image/', fileAbr ...
        ,' registered.tif');

    img = imread(imgFilename);
    % Note that the img pixels is Y x X, so the intensity is img(y, x)

    imshow(2 * mat2gray(img));
    hold on;

    % 2. Read the x, y coordinates of 81 reference dots 
    cordFilename = ['/Users/jackie/Downloads/SGTC microscope/' ...
        '2016-07-22-81 strands Atto sequence color/' ...
        'Analyze/axisAmineG25G2s.csv'];
    coord = csvread(cordFilename, 1, 0);

    % refSpot1 = refCoord(spot1Index,2:3);
    % refSpot2 = refCoord(spot2Index,2:3);
    % coord = regidRotate(refCoord(:, 2:3), refSpot1, refSpot2, spot1, spot2);
    % coord = [refCoord(:,1) coord];

    radius = 18;
    bgMinRadius = 25;
    bgMaxRadius = 28;
    numSpots = size(coord,1);
    
    for i = 1 : numSpots
        x0 = coord(i,2);
        y0 = coord(i,3);
        plot(x0, y0, 'r', 'MarkerSize', 100);
        coord(i,4) = averageIntensity(img, x0, y0, 0, radius);
        coord(i,5) = averageIntensity(img, x0, y0, bgMinRadius, bgMaxRadius);
    end

    intensityFilename = strcat('/Users/jackie/Downloads/SGTC microscope/' ...
        ,'2016-07-22-81 strands Atto sequence color/' ...
        ,'Analyze/intensity ', fileAbr, '.csv');    

    fid = fopen(intensityFilename,'w');
    fprintf(fid,'index,xcoord,ycoord,avgIntensity,bgIntensity\n');
    fclose(fid);
    dlmwrite(intensityFilename,coord,'-append','delimiter',',');
    %csvwrite(intensityFilename, coord);
end






