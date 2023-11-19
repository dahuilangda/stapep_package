from transformers import AutoTokenizer, EsmForProteinFolding
import torch


def predict_pdb(seq, device="cpu"):
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)

    device = torch.device(device)
    model.esm = model.esm.float()
    model = model.to(device)
    model.trunk.set_chunk_size(64)


    tokenized_input = tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"]
    tokenized_input = tokenized_input.to(device)

    with torch.no_grad():
        prediction = model.infer_pdb(seq)
    # with open(output, "w") as f:
    #     f.write(prediction)
    return prediction


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-seq", help="amino acid sequence", required=True)
    parser.add_argument("-device", help="device", default="cpu")
    args = parser.parse_args()

    print(predict_pdb(args.seq, args.device))
    # python -m stapep.esmfold -seq "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG" > test.pdb