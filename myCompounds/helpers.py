from compounddb.models import Compound, Tag
from targetsearch.helpers import get_chembl_sdfs, get_chembl_ids_with_struct, get_chembl_smiles
import time
from tools.models import Job
from tools.runapp import createJob, updateJob
from base64 import b64encode
from pubchem_rest_interface.Pubchem_pug import pubchemDownload

def addCompoundsAjax(user, source_id, ids, tags):
    if len(tags) != 0:
        Tag.ensureAllExist(tags, user)

    if len(ids) == 0:
        raise Exception('Empty list of "ids".')

    if source_id == 'chembl':
        sdfs = get_chembl_sdfs(ids)
        sdf = '\n$$$$\n'.join(sdfs) + '\n$$$$\n'

        newJob = createJob(user, 'Upload Compounds', '',
                            ['--user=' + str(user.id), '--tags=' + (','.join(tags))], sdf)
        time.sleep(2)
        newJob = updateJob(user, newJob.id)
        if newJob.status == Job.RUNNING:
            ajaxResponse = { 'success': True, 'message': 'Compound upload in progress. Check "Past Jobs" for status.' }
        elif newJob.status == Job.FINISHED:
            ajaxResponse = { 'success': True, 'message': 'Compounds uploaded successfully.'}
        else:
            ajaxResponse = { 'success': False, 'message': 'An error occurred while uploading your compounds.' }
        return ajaxResponse
    else:
        raise Exception('Unknown source_id: {}'.format(source_id))

def checkCompoundsAjax(user, source_id, ids):
    """Pre-checks a list of IDs by sorting them into lists of acceptable and
    problematic IDs. Problematic compounds include those that already exist in
    the workbench, and those that do not have structure information available
    in the ChEMBL database."""

    if len(ids) == 0:
        raise Exception('Empty list of "ids".')

    results = list()

    results_accept = dict()
    results_accept['select'] = True
    results_accept['enable'] = True
    results_accept['desc'] = "The following compounds can be added to your workbench."
    results_accept['ids'] = list()

    results_inWB = dict()
    results_inWB['select'] = False
    results_inWB['enable'] = True
    results_inWB['desc'] = "The following compounds are already in your workbench. You may still import them by checking the box next to the compound. (Only the ChEMBL ID was checked for duplicates.)"
    results_inWB['ids'] = list()

    results_nostruct = dict()
    results_nostruct['select'] = False
    results_nostruct['enable'] = False
    results_nostruct['desc'] = "The following compounds can't be imported to the workbench because compound structure information is missing from the current ChEMBL database. If you have this information, you may import these compounds manually (i.e. via SDF upload)."
    results_nostruct['ids'] = list()

    if source_id == 'chembl':
        ids_with_struct = get_chembl_ids_with_struct(ids)

        results_nostruct['ids'] = list(set(ids) - set(ids_with_struct))

        inWB = Compound.inWorkbench(user, ids_with_struct)
        for r in ids_with_struct:
            if r in inWB:
                results_inWB['ids'].append(r)
            else:
                results_accept['ids'].append(r)

        if len(results_accept['ids']) != 0:
            results.append(results_accept)
        if len(results_inWB['ids']) != 0:
            results.append(results_inWB)
        if len(results_nostruct['ids']) != 0:
            results.append(results_nostruct)

        return { 'success': True, 'results': results }
    else:
        raise Exception('Unknown source_id: {}'.format(source_id))

def downloadCompoundsAjax(user, source_id, ids, output_format, tags):
    """AJAX version of downloadCompounds(). Given a list of compound IDs
    (i.e. ChEMBL, PubChem, workbench ID, tags...), download the associated
    compound data (i.e. SDF, SMILES). The data will be returned in Base64."""

    if len(ids) == 0:
        raise Exception('Empty list of "ids".')

    if source_id == "chembl":
        if output_format == "sdf":
            sdfs = get_chembl_sdfs(tuple(ids))
            sdf = "\n$$$$\n".join(sdfs) + "\n$$$$\n"
            return { "success": True, "data": b64encode(sdf.encode()) }
        elif output_format == "smi":
            smiles = get_chembl_smiles(tuple(ids))
            smi = "\n".join(smiles) + "\n"
            return { "success": True, "data": b64encode(smi.encode()) }
        else:
            raise Exception("Invalid output_format: {}".format(output_format))
    elif source_id == "cid":
        raise Exception("Implement me...")
    elif source_id == "tag":
        raise Exception("Implement me...")
    elif source_id == "pubchem":
        raise Exception("Implement me...")
    else:
        raise Exception("Unknown source_id: {}".format(source_id))
