from pycisTopic.lda_models import *
import os
import sys
os.environ["MALLET_MEMORY"] = "100G"

def main(n_topic, wdir):
    cistopic_obj = pickle.load(open(f"{wdir}/pycistopic_outputs/cistopic_object.pkl", "rb"))
    out_dir = f"{wdir}/pycistopic_outputs/model_out/"
    tmp_path = "{wdir}/.tmp/"
    path_to_mallet_binary = "/tscc/nfs/home/rlancione/software/Mallet-202108/bin/mallet"
    os.makedirs(tmp_path, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    model = run_cgs_models_mallet(
        mallet_path=path_to_mallet_binary,
        cistopic_obj=cistopic_obj,
        n_topics=[int(n_topic)],
        n_cpu=14,
        n_iter=500,
        random_state=555,
        alpha=50,
        alpha_by_topic=True,
        eta=0.1,
        eta_by_topic=False,
        tmp_path=tmp_path,
        save_path=out_dir
    )

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
