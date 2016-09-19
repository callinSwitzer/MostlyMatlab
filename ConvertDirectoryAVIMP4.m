%% Make movie into mp4 instead of avi to put into keynote
%% choose a file in the directory you want
[FileName,PathName,FilterIndex] = uigetfile({'*.avi'},'File Selector');

%% Get vidnames
s = dir(fullfile(PathName,'*.avi'));
s = s(arrayfun(@(x) ~strcmp(x.name(1),'.'),s)) % ignore invisible files
vidNames = {s.name}';

%%

for ii = 4%length(vidNames)
    FileName = vidNames(ii); 
    % where you want to save it
    %file = '/Users/callinswitzer/Desktop/btag0.mp4' 
    foo = strsplit(char(FileName), '.')

    file = char(strcat(PathName, foo(1), '.mp4'));

    vidObj = VideoReader(strcat(PathName, char(FileName)));

    % get info
    vidHeight = vidObj.Height;
    vidWidth = vidObj.Width;

    s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
        'colormap',[]);

    writer = VideoWriter(file, 'MPEG-4');

    writer.FrameRate = vidObj.FrameRate;

    open(writer);
    % Read and write each frame.


    % read in frames
    k = 1;
    while hasFrame(vidObj)
        s(k).cdata = readFrame(vidObj);
        
        % watermark
        positions = [round(3.5/5 * size(s(k).cdata, 1)) round(3.5/5 * size(s(k).cdata, 2)); round(3.5/5 * size(s(k).cdata, 1)) round(3.5/5 * size(s(k).cdata, 2) + 40)];
        test_str = cell({'Randolph VT' ,'Callin Switzer'});

        RGB = insertText(s(k).cdata, positions,...
            test_str,'BoxOpacity',0, 'TextColor', [150, 150, 150],...
            'FontSize', 30);
        %imshow(imadjust(RGB,[;],[], .6))
   
        writeVideo(writer,imadjust(RGB,[;],[], .6));
        k = k+1
    end
    close(writer);
end