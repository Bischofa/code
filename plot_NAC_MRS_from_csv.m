%%
% included voxels
fname = '/data/waubant1/7T_NAC_MS_GSH/3DMRSI_results/7T_NAC_GSH_MRSI_AtlasROI_metabolite_ratios_cut20_20180723.csv';
% excluded voxels
% fname = '7T_NAC_GSH_MRSI_AtlasROI_metabolites_cut20_badvoxellist_20180723_excluded_voxels.csv'

addpath('/netopt/share/lib/local/brain/matlab/');

metmat = readtable(fname);

dataroot = '/data/waubant1/7T_NAC_MS_GSH';
lcquantroot = 'comb_cor_sum';
exam_pre = {' '};

n = 1;
while n <= size(metmat,1)

    bnum = metmat{n,2};
    tnum = metmat{n,3};
    exam = metmat{n,4};
    column = metmat{n,8};
    row = metmat{n,9};
    slice = metmat{n,10};
    rvalue = metmat{n,5};
    hemi = metmat{n,6};
    roiloc = [column row slice];
    
    dstpath = sprintf('%s/%s/%s',dataroot,bnum{:},tnum{:});
    if strcmp(exam_pre{:},exam{:}) == 0
        t1vimage = read_idf_image(sprintf('%s/spectra_csi/%s_t1v_resampled',dstpath,exam{:}));
        temproiname = sprintf('%s/HOAtlas_analysis/HarvardOxford_analysis/ROI_Spec/HOcort25_to_%s_t1v_val29_lh_spec_percent',...
            dstpath,exam{:});
        roitemp = read_idf_image(temproiname);
        speccor = read_ddf_image(sprintf('%s/spectra_csi/LcGrid/%s_csi_%s_LCM_cor',dstpath,exam{:},lcquantroot));
        specbaseline = read_ddf_image(sprintf('%s/spectra_csi/LcGrid/%s_csi_%s_LCM_baseline',dstpath,exam{:},lcquantroot));
        specphase = read_ddf_image(sprintf('%s/spectra_csi/LcGrid/%s_csi_%s_LCM_phased',dstpath,exam{:},lcquantroot));
    end
    
    plotzoom = (speccor.ddf.ppm_ref-[4.1 1.8]).*speccor.ddf.centfreq;
    plotzoom = plotzoom*(speccor.ddf.specpoints/speccor.ddf.sweepwidth) + speccor.ddf.specpoints/2;
    plotzoom = round(plotzoom);
    figure(1);
    subplot(1,3,1);
    plot(real(speccor.img(plotzoom(1):plotzoom(2),roiloc(1),roiloc(2),roiloc(3))));hold on
    plot(real(specbaseline.img(plotzoom(1):plotzoom(2),roiloc(1),roiloc(2),roiloc(3))));
    plot(real(specphase.img(plotzoom(1):plotzoom(2),roiloc(1),roiloc(2),roiloc(3))));hold off
    title(sprintf('[%d %d %d] %d/%d)',roiloc(1),roiloc(2),roiloc(3),n,size(metmat,1)));
    
    
    
    subplot(1,3,3);
    roispace = zeros(3,1);
    roispace(:) = NaN;
    for i = 1:3
        a = roitemp.idf.dcos(i,:).*((roiloc-1).*roitemp.idf.pixelsize);
        roispace(i) = roitemp.idf.toplc(i)+a(1)+a(2)+a(3);
    end
    B = roispace-t1vimage.idf.toplc';
    A = [t1vimage.idf.dcos(1,:).*t1vimage.idf.pixelsize;t1vimage.idf.dcos(2,:).*t1vimage.idf.pixelsize;t1vimage.idf.dcos(3,:).*t1vimage.idf.pixelsize];
    roit1vloc = floor(A\B);
    
    roiimage = ones(size(t1vimage.img,1),size(t1vimage.img,2));
    a= floor([roitemp.idf.pixelsize(1) roitemp.idf.pixelsize(2)]./t1vimage.idf.pixelsize(1:2));
    roiimage(roit1vloc(1)-a(1)/2:roit1vloc(1)+a(1)/2,roit1vloc(2)-a(2)/2:roit1vloc(2)+a(2)/2) =0;
    imageslice = t1vimage.img(:,:,roit1vloc(3,1))';
    h=imshow(imageslice,[min(imageslice(:)) max(imageslice(:))]);
    set(h, 'AlphaData', roiimage');
    title(sprintf('%s %s %s %s',bnum{:},exam{:},rvalue{:},hemi{:}));
    
    subplot(1,3,2);
    imshow(imageslice,[min(imageslice(:)) max(imageslice(:))]);
    
   
    prompt = 'Next Display? negative = previous  (e.g. -1 = previous voxel); enter = next voxel; positive = next (e.g. +10 = 10 voxels later)  ';
    voxelloc = input(prompt);
    if length(voxelloc)> 0       
        n = n + voxelloc;
    else
        n = n+1;
    end
    
    exam_pre = exam;

end



