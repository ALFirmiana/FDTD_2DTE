import h5py
import matplotlib.pyplot as plt

file_path = "build/test.h5"

# |%%--%%| <aUoEO2UEoD|1hQv8kKawq>


def print_h5_items(name, obj):
    print(name)
    if isinstance(obj, h5py.Dataset):
        print("    Type: Dataset")
        print("    Shape:", obj.shape)
        print("    Data type:", obj.dtype)
    elif isinstance(obj, h5py.Group):
        print("    Type: Group")

    # 打印属性
    for key, value in obj.attrs.items():
        print("    Attribute:", key, "=", value)


# 打开HDF5文件
with h5py.File(file_path, 'r') as file:
    # 遍历文件内容
    file.visititems(print_h5_items)
# |%%--%%| <1hQv8kKawq|OjNewbU3fY>

with h5py.File(file_path, 'r') as file:
    plt.imshow(file["Bz_9"])

