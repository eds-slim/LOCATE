function copying_image_geometry(source_file, target_file, verbose)

setenv( 'FSLDIR', '/usr/local/fsl' );
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

setenv('PATH', [getenv('PATH'),':/usr/local/fsl/bin'])

[status,output] = call_fsl(sprintf('fslcpgeom %s %s',source_file, target_file));

if verbose
    if ~status
        fprintf('Successfully copied the image dimensions... \n');
    else
        fprintf(sprintf('Failed copying dimensions from %s to %s. Kindly check your FSL installation. If everything is ok, try copying dimesnions using fslcpgeom',...
            source_file,target_file));
    end
end
        

