import os, glob, shutil, sys, inspect, errno

from distutils.spawn import find_executable
from mrtrix3 import app, fsl, file, image, path, run


app.init ('Lea Vinokur (lea.vinokur@gmail.com)',
          'Perform group analysis of diffusion MRI data with a Fixel-Based Analysis (FBA) of Fibre Density, Fibre Cross-section and a combined measure (Fibre Density & Cross-section). The analysis pipeline relies primarily on the MRtrix3 software package (www.mrtrix.org).')

app.cmdline.add_argument('in_dir', help='The directory with the input dataset '
                                        'formatted according to the BIDS standard.')
app.cmdline.add_argument('output_dir', help='The directory where the output files '
                                         'should be stored. If you are running group level analysis '
                                         'this folder should be prepopulated with the results of the '
                                         'participant level analysis.')

analysis_level_choices = ['participant1', 'group1', 'participant2', 'group2']

app.cmdline.add_argument('analysis_level', help='Level of the analysis that will be performed. '
                                                   'Valid choices are: [' + ', '.join(analysis_level_choices) + ']. \nMultiple participant '
                                                   'level analyses can be run independently(in parallel) using the same output_dir.',
                                              choices = analysis_level_choices)

app.cmdline.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label '
                                                 'corresponds to sub-<participant_label> from the BIDS spec '
                                                 '(so it does not include "sub-"). If this parameter is not '
                                                 'provided all subjects should be analyzed. Multiple '
                                                 'participants can be specified with a space separated list.',
                                                  nargs='+')
app.parse()

if app.isWindows():
    app.error('Script cannot be run on Windows due to FSL dependency')

subjects_to_analyze = []
# only for a subset of subjects
if app.args.participant_label:
    subjects_to_analyze = app.args.participant_label
# for all subjects
else:
    subject_dirs = glob.glob(os.path.join(app.args.in_dir, 'sub-*'))
    subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

# create output subjects directory
all_subjects_dir = os.path.join(app.args.output_dir, 'subjects');
if not os.path.exists(all_subjects_dir):
    os.mkdir(all_subjects_dir)

# create the output template directory
template_dir = os.path.join(app.args.output_dir, 'template')
if not os.path.exists(template_dir):
    os.mkdir(template_dir)

# create a temporary directory for intermediate files
app.makeTempDir()



# loop over subjects and perform subject-level analysis (genereta mif, biascorrect, response functions)
if app.args.analysis_level == 'participant1':
    for subject_label in subjects_to_analyze:
        label = 'sub-' + subject_label
        print('running basic pre-processing for ' + label)

        # Read DWI(s) in BIDS folder
        all_dwi_images = glob.glob(os.path.join(app.args.in_dir, label, '*dwi', '*_dwi.nii*'))

        if (len(all_dwi_images) > 1):
            app.error('Detecting more that one DWI - please fix')

        # Create output subject directory
        subject_dir = os.path.join(all_subjects_dir, subject_label)
        if not os.path.exists(subject_dir):
            os.mkdir(subject_dir)

        # Check existence output files from this analysis level
        dwi_preproc_file = os.path.join(subject_dir, 'dwi_preproc.mif')
        app.checkOutputPath(dwi_preproc_file)
        wm_response_file = os.path.join(subject_dir, 'wm_response.txt')
        app.checkOutputPath(wm_response_file)
        gm_response_file = os.path.join(subject_dir, 'gm_response.txt')
        app.checkOutputPath(gm_response_file)
        csf_response_file = os.path.join(subject_dir, 'csf_response.txt')
        app.checkOutputPath(csf_response_file)

        # DW gradient files
        grad_prefix = os.path.join(app.args.in_dir, label, 'dwi', label + '_dwi')
        if not (os.path.isfile(grad_prefix + '.bval') and os.path.isfile(grad_prefix + '.bvec')):
            grad_prefix = os.path.join(app.args.in_dir, 'dwi')
            if not (os.path.isfile(grad_prefix + '.bval') and os.path.isfile(grad_prefix + '.bvec')):
                app.error('Unable to locate valid diffusion gradient table');
        grad_import_option = ' -fslgrad ' + grad_prefix + '.bvec ' + grad_prefix + '.bval'

        # Stuff DWI gradients in *.mif file
        dwi_mrtrix_file = os.path.join(app.tempDir, subject_label + '_dwi.mif')
        print(dwi_mrtrix_file)
        # run.command('mrconvert ' + all_dwi_images[0] + grad_import_option + ' ' + dwi_mrtrix_file)
        # dwi_biascorrect_file = os.path.join(app.tempDir, subject_label + 'dwi_biascorr.mif')
        print('mrconvert command ...... ')
        print('mrconvert ' + all_dwi_images[0] + grad_import_option + ' ' + dwi_mrtrix_file)
        dwi_biascorrect_file = os.path.join(app.tempDir, subject_label + '_dwi_biascorr.mif')
        # run.command('dwibiascorrect' + ' -ants ' + dwi_mrtrix_file + ' ' + dwi_denoised_file)
        print('bias correction ......... ')
        print('dwibiascorrect' + ' -ants ' + dwi_mrtrix_file + ' ' + dwi_biascorrect_file)
        # Estimate WM, GM, CSF response functions
        # run.command('dwi2response ' + app.numThreads + ' dhollander ' + dwi_preproc_file + ' ' + wm_response_file + ' '
        #                                                                                           + gm_response_file + ' '
        #                                                                                           + csf_response_file + app.force)
        print('response estimation ..... ')
        print('dwi2response ' + ' dhollander ' + dwi_biascorrect_file + ' ' + wm_response_file + ' ' + gm_response_file + ' ' + csf_response_file)

# running group level 1 (average response functions) TODO check for user supplied subset to ensure response is not biased
elif app.args.analysis_level == "group1":
    print('averaging response functions')

    # Check output files exist
    wm_response_file = os.path.join(app.args.output_dir, 'average_wm_response.txt')
    app.checkOutputPath(wm_response_file)
    gm_response_file = os.path.join(app.args.output_dir, 'average_gm_response.txt')
    app.checkOutputPath(gm_response_file)
    csf_response_file = os.path.join(app.args.output_dir, 'average_csf_response.txt')
    app.checkOutputPath(csf_response_file)

    # process subset
    input_wm_files = []
    input_gm_files = []
    input_csf_files = []

    input_wm_files = glob.glob(os.path.join(all_subjects_dir, '*', 'wm_response.txt'))
    input_gm_files = glob.glob(os.path.join(all_subjects_dir, '*', 'gm_response.txt'))
    input_csf_files = glob.glob(os.path.join(all_subjects_dir, '*', 'csf_response.txt'))

    run.command('average_response ' + ' '.join(input_wm_files) + ' ' + wm_response_file + app.force)
    run.command('average_response ' + ' '.join(input_gm_files) + ' ' + gm_response_file + app.force)
    run.command('average_response ' + ' '.join(input_csf_files) + ' ' + csf_response_file + app.force)


# running participant level 2 (upsample, compute brain masks and FODs, perform intensity normalisation and bias field correction)
elif app.args.analysis_level == "participant2":
    print ('calculating masks and performing CSD')
    for subject_label in subjects_to_analyze:

        subject_dir = os.path.join(all_subjects_dir, subject_label)
        output_mask = os.path.join(subject_dir, 'mask.mif')
        app.checkOutputPath(output_mask)
        output_fod = os.path.join(subject_dir, 'fod.mif')
        app.checkOutputPath(output_fod)
        output_gm = os.path.join(subject_dir, 'gm.mif')
        app.checkOutputPath(output_gm)
        output_csf = os.path.join(subject_dir, 'csf.mif')
        app.checkOutputPath(output_csf)

        min_voxel_size = 1.25;
        if app.args.vox_size:
          min_voxel_size = float(app.args.vox_size)

        voxel_sizes = getHeaderInfo(os.path.join(subject_dir, 'dwi_preproc.mif'), 'vox').split()
        mean_voxel_size = 0.0
        for i in range(0,3):
          mean_voxel_size = mean_voxel_size + float(voxel_sizes[i]) / 3.0

        input_to_csd = os.path.join(subject_dir, 'dwi_preproc.mif')
    # Compute brain mask
    run.command('dwi2mask ' + input_to_csd + ' ' + output_mask + app.force)

    #TODO add dilate mask by a voxel or two

    # Perform CSD
    run.command('dwi2fod msmt_csd ' + input_to_csd + ' -mask ' + output_mask + ' ' +
                os.path.join(app.args.output_dir, 'average_wm_response.txt') + ' ' +  os.path.join(app.tempDir, subject_label + 'fod.mif') + ' ' +
                os.path.join(app.args.output_dir, 'average_gm_response.txt') + ' ' + os.path.join(app.tempDir, subject_label + 'gm.mif') + ' ' +
                os.path.join(app.args.output_dir, 'average_csf_response.txt') + ' ' + os.path.join(app.tempDir, subject_label + 'csf.mif'))

    run.command('mtbin -independent ' + os.path.join(app.tempDir, subject_label + 'fod.mif') + ' ' + output_fod + ' ' +
                                       os.path.join(app.tempDir, subject_label + 'gm.mif')  + ' ' + output_gm  + ' ' +
                                       os.path.join(app.tempDir, subject_label + 'csf.mif') + ' ' + output_csf + app.force)


# running group level 2 (generate FOD template)
elif app.args.analysis_level == 'group2':

    # Check if outputs exist
    fod_template = os.path.join(template_dir, 'fod_template.mif')
    app.checkOutputPath(fod_template)

    app.gotoTempDir()
    os.mkdir('fod_input')
    os.mkdir('mask_input')

    # Check if all members of subset exist
    if (len(subset) > 0):
        print('Using a group subset to compute population template' + str(subset))
        subject_labels = [os.path.basename(x) for x in glob.glob(os.path.join(all_subjects_dir, '*'))]
    for subj in subset:
        if subj not in subject_labels:
            app.error('subject label (' + os.path.basename(subj) + ') supplied as part of -group_subset option does exist in subjects directory')
            app.checkOutputPath(os.path.join(all_subjects_dir, subj, 'subject2template_warp.mif'))
            app.checkOutputPath(os.path.join(all_subjects_dir, subj, 'template2subject_warp.mif'))


    # make symlinks to all population_template inputs in single directory
    for subj in glob.glob(os.path.join(all_subjects_dir, '*')):
        if (len(subset) > 0):
            if os.path.basename(subj) in subset:
                os.symlink(os.path.join(subj, 'fod.mif'), os.path.join('fod_input', os.path.basename(subj) + '.mif'))
                os.symlink(os.path.join(subj, 'mask.mif'), os.path.join('mask_input', os.path.basename(subj) + '.mif'))
        else:
          os.symlink(os.path.join(subj, 'fod.mif'), os.path.join('fod_input', os.path.basename(subj) + '.mif'))
          os.symlink(os.path.join(subj, 'mask.mif'), os.path.join('mask_input', os.path.basename(subj) + '.mif'))

    # Compute FOD template
    if (len(subset) > 0):
        run.command('population_template fod_input -mask mask_input ' + os.path.join(app.tempDir, 'tmp.mif'))
    # Set a field in the header of the template to mark it as being generated as a subset or not. This is used in the next step
        run.command('mrconvert ' + os.path.join(app.tempDir, 'tmp.mif') + ' -header_set made_from_subset true ' + fod_template + app.force)
    else:
        run.command('population_template fod_input -mask mask_input ' + os.path.join(app.tempDir, 'tmp.mif') + ' -warp_dir ' +  os.path.join(app.tempDir, 'warps'))
        run.command('mrconvert ' + os.path.join(app.tempDir, 'tmp.mif') + ' -header_set made_from_subset false ' + fod_template + app.force)
    # Save all warps since we don't need to generate them in the next step if all subjects were used to make the template
    for subj in [os.path.basename(x) for x in glob.glob(os.path.join(all_subjects_dir, '*'))]:
        run.command('warpconvert -type warpfull2deformation -template ' + fod_template + ' '
                              + os.path.join(app.tempDir, 'warps', subj + '.mif') + ' '
                              + os.path.join(all_subjects_dir, subj, 'subject2template_warp.mif') + app.force)
        run.command('warpconvert -type warpfull2deformation -from 2 -template ' + os.path.join(all_subjects_dir, subj, 'fod.mif') + ' '
                              + os.path.join(app.tempDir, 'warps', subj + '.mif') + ' '
                              + os.path.join(all_subjects_dir, subj, 'template2subject_warp.mif') + app.force)
