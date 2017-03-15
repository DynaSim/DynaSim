function varargout = mkdir_silent (output_path,varargin)
    % Deletes folder contents then recreates

%     fprintf(['Deleting all data in folder: \n' output_path '.\nPress any key to continue or CTRL-C to abort.\n']);
%     pause

    suppress_output = 1;

    if ~exist(output_path,'file')
        fprintf('Creating %s \n', output_path);
        %system( ['mkdir ' output_path]);
        [varargout{1:nargout}] = mkdir(output_path,varargin{:});   % Outputs are supplied simply to suppress warning
                                        % message for existing folder.
    else
        if ~suppress_output; fprintf('Folder %s already exists. Doing nothing.\n',output_path); end
    end
    
end