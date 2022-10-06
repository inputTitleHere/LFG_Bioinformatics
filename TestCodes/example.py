print("Hello world!")
print("Print this for me!")



cancer=['Brca','LUAD',"USCS","BLCA"]
dic = {}
target_gene=["CTGF","CYR61"]

for i in cancer:
    print(f"Current cancer name is {i}")
    print("Processing somthibng")
    for k in target_gene:
        print(f"doing {i} on {k}")
