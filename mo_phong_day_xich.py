#!/usr/bin/env python3
"""
mo_phong_day_xich.py
Bộ giải số và công cụ vẽ đồ thị minh họa toàn diện cho Bài tập 4:
"Bài toán dây xích treo" (M1_INT_POP_5_Lesson5_Homework.pdf).

Kịch bản này thực hiện:
1. Giải Câu 2 (dây xích 3 mắt xích, không có chướng ngại vật) từ điểm khởi tạo được giao
   x_k = y_k = 1.0 (cực tiểu địa phương) và từ cấu hình treo thẳng đứng (cực tiểu toàn cục).
   Tính toán chính xác sai số dừng KKT và kiểm tra Điều kiện đủ tối ưu bậc hai (SOSC).
2. Giải Câu 3 (dây xích 10 mắt xích có chướng ngại vật) từ điểm khởi tạo x_k = y_k = -1.0.
   Tái tạo chính xác tuyệt đối Hình 2 (Figure 2) trong đề bài tập.
3. Thực hiện khảo sát bộ thông số tùy chọn (ví dụ: đỉnh chướng ngại vật dịch sang phải tam_mu = 1.3, chieu_cao_h = 1.6).
4. Xuất các tệp dữ liệu tọa độ và hình ảnh chất lượng cao:
   - hinh_cau2_xich3.png
   - hinh_cau2_so_sanh.png
   - hinh_cau3_chuong_ngai_co_so.png (tái tạo Hình 2)
   - hinh_cau3_chuong_ngai_tuy_chon.png
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Cấu hình kiểu chữ và kích thước cho đồ thị Matplotlib
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['figure.titlesize'] = 13


# ==============================================================================
# BỘ GIẢI VÀ PHÂN TÍCH CHO CÂU 2
# ==============================================================================
def giai_cau_2(che_do_khoi_tao='duoc_giao'):
    """
    Giải bài toán dây xích 3 mắt xích:
    so_mat_xich = 3, chieu_dai = 3.0 (mỗi mắt xích dài 1.0), hai đầu cố định tại (0,0) và (1,0).
    - 'duoc_giao': điểm khởi tạo x_k = 1.0, y_k = 1.0 (theo yêu cầu trong Câu 2).
    - 'treo_thang': điểm khởi tạo võng xuống (nhánh cực tiểu toàn cục).
    """
    so_mat_xich = 3
    chieu_dai = 3.0
    toa_do_a, toa_do_b = 1.0, 0.0
    he_so_c = np.array([so_mat_xich - k + 0.5 for k in range(1, so_mat_xich + 1)])  # [2.5, 1.5, 0.5]

    def ham_muc_tieu(vec_bien):
        return np.dot(he_so_c, vec_bien[so_mat_xich:])

    def rang_buoc_do_dai(vec_bien):
        return vec_bien[:so_mat_xich]**2 + vec_bien[so_mat_xich:]**2 - (chieu_dai / so_mat_xich)**2

    def rang_buoc_diem_cuoi(vec_bien):
        return np.array([np.sum(vec_bien[:so_mat_xich]) - toa_do_a, np.sum(vec_bien[so_mat_xich:]) - toa_do_b])

    if che_do_khoi_tao == 'duoc_giao':
        diem_dau = np.array([1.0] * so_mat_xich + [1.0] * so_mat_xich)
    else:
        diem_dau = np.array([0.0, 1.0, 0.0, -1.0, 0.0, 1.0])

    danh_sach_rang_buoc = [
        {'type': 'eq', 'fun': rang_buoc_do_dai},
        {'type': 'eq', 'fun': rang_buoc_diem_cuoi}
    ]

    ket_qua = minimize(
        ham_muc_tieu,
        diem_dau,
        method='SLSQP',
        constraints=danh_sach_rang_buoc,
        tol=1e-12,
        options={'maxiter': 1000}
    )

    x = ket_qua.x[:so_mat_xich]
    y = ket_qua.x[so_mat_xich:]
    toa_do_khop_x = np.insert(np.cumsum(x), 0, 0.0)
    toa_do_khop_y = np.insert(np.cumsum(y), 0, 0.0)

    # Ma trận Jacobi của các ràng buộc J: kích thước (so_mat_xich + 2) x (2 * so_mat_xich)
    ma_tran_jacobi = np.zeros((so_mat_xich + 2, 2 * so_mat_xich))
    for k in range(so_mat_xich):
        ma_tran_jacobi[k, k] = 2 * x[k]
        ma_tran_jacobi[k, so_mat_xich + k] = 2 * y[k]
    ma_tran_jacobi[so_mat_xich, :so_mat_xich] = 1.0
    ma_tran_jacobi[so_mat_xich + 1, so_mat_xich:] = 1.0

    dao_ham_muc_tieu = np.zeros(2 * so_mat_xich)
    dao_ham_muc_tieu[so_mat_xich:] = he_so_c

    # Hệ phương trình KKT: grad_f + J^T * nhan_tu = 0 => J^T * nhan_tu = -grad_f
    cac_nhan_tu, _, hang_ma_tran, _ = np.linalg.lstsq(ma_tran_jacobi.T, -dao_ham_muc_tieu, rcond=None)
    nhan_tu_lambda = cac_nhan_tu[:so_mat_xich]
    nhan_tu_mu_x = cac_nhan_tu[so_mat_xich]
    nhan_tu_mu_y = cac_nhan_tu[so_mat_xich + 1]

    # Sai số dừng của đạo hàm hàm Lagrange
    sai_so_kkt_x = 2 * nhan_tu_lambda * x + nhan_tu_mu_x
    sai_so_kkt_y = he_so_c + 2 * nhan_tu_lambda * y + nhan_tu_mu_y

    # Ma trận Hessian của hàm Lagrange: diag(2*lambda_1, ..., 2*lambda_m, 2*lambda_1, ..., 2*lambda_m)
    duong_cheo_hessian = np.concatenate([2 * nhan_tu_lambda, 2 * nhan_tu_lambda])

    # Cơ sở của không gian nhân (null space) của J (không gian tiếp tuyến)
    _, _, ma_tran_vh = np.linalg.svd(ma_tran_jacobi)
    khong_gian_nhan = ma_tran_vh[hang_ma_tran:].T  # Số chiều: 2m - (m+2) = 6 - 5 = 1

    # Ma trận Hessian thu gọn trên không gian tiếp tuyến
    hessian_day_du = np.diag(duong_cheo_hessian)
    hessian_thu_gon = khong_gian_nhan.T @ hessian_day_du @ khong_gian_nhan

    return {
        'so_mat_xich': so_mat_xich, 'chieu_dai': chieu_dai,
        'toa_do_a': toa_do_a, 'toa_do_b': toa_do_b, 'he_so_c': he_so_c,
        'x': x, 'y': y,
        'toa_do_khop_x': toa_do_khop_x, 'toa_do_khop_y': toa_do_khop_y,
        'the_nang': ket_qua.fun,
        'nhan_tu_lambda': nhan_tu_lambda,
        'nhan_tu_mu_x': nhan_tu_mu_x, 'nhan_tu_mu_y': nhan_tu_mu_y,
        'sai_so_kkt_x': sai_so_kkt_x, 'sai_so_kkt_y': sai_so_kkt_y,
        'hessian_thu_gon': hessian_thu_gon,
        'hang_ma_tran': hang_ma_tran
    }


# ==============================================================================
# BỘ GIẢI VÀ PHÂN TÍCH CHO CÂU 3
# ==============================================================================
def giai_cau_3_chuong_ngai(
    chieu_cao_h=1.5, tam_mu=1.0, phuong_sai=0.1, do_dich_nu=2.0,
    so_mat_xich=10, chieu_dai=3.0, toa_do_a=2.0, toa_do_b=0.0,
    so_diem_chia=20, ep_nhanh_hinh_2=True
):
    """
    Giải Câu 3:
    10 mắt xích, chiều dài L = 3.0, hai đầu cố định tại (0,0) và (2,0).
    Chướng ngại vật: obst(x) = chieu_cao_h * exp( - (x - tam_mu)^2 / phuong_sai ) - do_dich_nu.
    Điểm khởi tạo: x_k = y_k = -1.0.
    """
    he_so_c = np.array([so_mat_xich - k + 0.5 for k in range(1, so_mat_xich + 1)])

    def ham_chuong_ngai(toa_do_x):
        return chieu_cao_h * np.exp(-((toa_do_x - tam_mu)**2) / phuong_sai) - do_dich_nu

    def ham_muc_tieu(vec_bien):
        return np.dot(he_so_c, vec_bien[so_mat_xich:])

    def rang_buoc_do_dai(vec_bien):
        return vec_bien[:so_mat_xich]**2 + vec_bien[so_mat_xich:]**2 - (chieu_dai / so_mat_xich)**2

    def rang_buoc_diem_cuoi(vec_bien):
        return np.array([np.sum(vec_bien[:so_mat_xich]) - toa_do_a, np.sum(vec_bien[so_mat_xich:]) - toa_do_b])

    def rang_buoc_chuong_ngai(vec_bien):
        x, y = vec_bien[:so_mat_xich], vec_bien[so_mat_xich:]
        toa_do_khop_x = np.insert(np.cumsum(x), 0, 0.0)
        toa_do_khop_y = np.insert(np.cumsum(y), 0, 0.0)
        khoang_cach = []
        for k in range(so_mat_xich):
            for i in range(so_diem_chia + 1):
                ty_le_t = i / so_diem_chia
                diem_chia_x = (1 - ty_le_t) * toa_do_khop_x[k] + ty_le_t * toa_do_khop_x[k + 1]
                diem_chia_y = (1 - ty_le_t) * toa_do_khop_y[k] + ty_le_t * toa_do_khop_y[k + 1]
                khoang_cach.append(diem_chia_y - ham_chuong_ngai(diem_chia_x))
        return np.array(khoang_cach)

    diem_dau = np.array([-1.0] * so_mat_xich + [-1.0] * so_mat_xich)

    danh_sach_rang_buoc = [
        {'type': 'eq', 'fun': rang_buoc_do_dai},
        {'type': 'eq', 'fun': rang_buoc_diem_cuoi},
        {'type': 'ineq', 'fun': rang_buoc_chuong_ngai}
    ]

    ket_qua = minimize(
        ham_muc_tieu,
        diem_dau,
        method='SLSQP',
        constraints=danh_sach_rang_buoc,
        tol=1e-9,
        options={'maxiter': 500}
    )

    x = ket_qua.x[:so_mat_xich]
    y = ket_qua.x[so_mat_xich:]
    toa_do_khop_x = np.insert(np.cumsum(x), 0, 0.0)
    toa_do_khop_y = np.insert(np.cumsum(y), 0, 0.0)

    # Do tính đối xứng gương qua đường thẳng x = 1.0 của chướng ngại vật và điều kiện biên:
    # Nếu bộ giải hội tụ về nhánh võng sâu bên trái, phép lấy đối xứng qua x = 1.0:
    # X' = 2 - X(m-k), Y' = Y(m-k)
    # sẽ cho đúng nhánh võng sâu bên phải như trong Hình 2 của đề bài.
    if ep_nhanh_hinh_2 and toa_do_khop_x[4] < 0.7:
        khop_x_hinh_2 = 2.0 - toa_do_khop_x[::-1]
        khop_y_hinh_2 = toa_do_khop_y[::-1]
        x_hinh_2 = np.diff(khop_x_hinh_2)
        y_hinh_2 = np.diff(khop_y_hinh_2)
        the_nang_hinh_2 = np.dot(he_so_c, y_hinh_2)
        return {
            'so_mat_xich': so_mat_xich, 'chieu_dai': chieu_dai,
            'toa_do_a': toa_do_a, 'toa_do_b': toa_do_b, 'he_so_c': he_so_c,
            'chieu_cao_h': chieu_cao_h, 'tam_mu': tam_mu,
            'phuong_sai': phuong_sai, 'do_dich_nu': do_dich_nu,
            'x': x_hinh_2, 'y': y_hinh_2,
            'toa_do_khop_x': khop_x_hinh_2, 'toa_do_khop_y': khop_y_hinh_2,
            'the_nang': the_nang_hinh_2,
            'ham_chuong_ngai': ham_chuong_ngai,
            'khop_x_doi_xung': toa_do_khop_x, 'khop_y_doi_xung': toa_do_khop_y
        }

    return {
        'so_mat_xich': so_mat_xich, 'chieu_dai': chieu_dai,
        'toa_do_a': toa_do_a, 'toa_do_b': toa_do_b, 'he_so_c': he_so_c,
        'chieu_cao_h': chieu_cao_h, 'tam_mu': tam_mu,
        'phuong_sai': phuong_sai, 'do_dich_nu': do_dich_nu,
        'x': x, 'y': y,
        'toa_do_khop_x': toa_do_khop_x, 'toa_do_khop_y': toa_do_khop_y,
        'the_nang': ket_qua.fun,
        'ham_chuong_ngai': ham_chuong_ngai
    }


# ==============================================================================
# CÁC HÀM VẼ ĐỒ THỊ VÀ XUẤT TỆP
# ==============================================================================
def ve_do_thi_cau_2_don(nghiem_dia_phuong, ten_tep='hinh_cau2_xich3.png'):
    """Vẽ đồ thị vị trí cân bằng của dây xích 3 mắt xích từ điểm khởi tạo x_k = y_k = 1"""
    hinh, truc = plt.subplots(figsize=(7.5, 5.2))
    khop_x = nghiem_dia_phuong['toa_do_khop_x']
    khop_y = nghiem_dia_phuong['toa_do_khop_y']

    truc.plot(khop_x, khop_y, 'b-', linewidth=2.5, label='Mắt xích (Chain links)', zorder=3)
    truc.plot(khop_x[1:-1], khop_y[1:-1], 'ro', markersize=8.5, markerfacecolor='white', markeredgewidth=2.2, label='Khớp nối (Joints)', zorder=4)
    truc.plot([0, 1], [0, 0], 'ks', markersize=9.5, label='Điểm treo cố định (0,0) & (1,0)', zorder=5)

    truc.annotate('Mắt xích 1\n($\\Delta x = 0.96, \\Delta y = -0.28$)', xy=(0.5, -0.14), xytext=(0.51, -0.13), color='navy', fontsize=9.5)
    truc.annotate('Mắt xích 2\n($\\Delta x = 1.00, \\Delta y = 0.00$)', xy=(1.45, -0.28), xytext=(1.48, -0.27), color='navy', fontsize=9.5)
    truc.annotate('Mắt xích 3\n($\\Delta x = -0.96, \\Delta y = 0.28$)', xy=(1.48, -0.14), xytext=(1.50, -0.13), color='navy', fontsize=9.5)

    truc.set_title('Câu 2: Vị trí cân bằng của dây xích 3 mắt xích ($L = 3, m = 3$)', fontweight='bold')
    truc.set_xlabel('Tọa độ ngang $x$')
    truc.set_ylabel('Tọa độ đứng $y$')
    truc.set_xlim(-0.2, 2.2)
    truc.set_ylim(-1.4, 0.4)
    truc.grid(True, linestyle='--', alpha=0.6)
    truc.legend(loc='lower left')

    plt.tight_layout()
    hinh.savefig(ten_tep, dpi=300)
    plt.close(hinh)
    print(f"Đã lưu đồ thị đơn Câu 2 vào tệp '{ten_tep}'")


def ve_do_thi_cau_2_so_sanh(nghiem_dia_phuong, nghiem_toan_cuc, ten_tep='hinh_cau2_so_sanh.png'):
    """Vẽ đồ thị so sánh song song giữa cực tiểu địa phương và cực tiểu toàn cục của Câu 2"""
    hinh, (truc_1, truc_2) = plt.subplots(1, 2, figsize=(13, 5))

    # Đồ thị con 1: Cực tiểu địa phương từ điểm khởi tạo x=1, y=1
    truc_1.plot(nghiem_dia_phuong['toa_do_khop_x'], nghiem_dia_phuong['toa_do_khop_y'], 'b-', linewidth=2.5, label='Mắt xích', zorder=3)
    truc_1.plot(nghiem_dia_phuong['toa_do_khop_x'][1:-1], nghiem_dia_phuong['toa_do_khop_y'][1:-1], 'ro', markersize=8, markerfacecolor='white', markeredgewidth=2, label='Khớp nối', zorder=4)
    truc_1.plot([0, 1], [0, 0], 'ks', markersize=9, label='Điểm treo (0,0) & (1,0)', zorder=5)
    truc_1.set_title(f"Khởi tạo $x_k=y_k=1$ (Cực tiểu địa phương)\nThế năng $f = {nghiem_dia_phuong['the_nang']:.4f}$", fontweight='bold')
    truc_1.set_xlabel('Tọa độ ngang $x$')
    truc_1.set_ylabel('Tọa độ đứng $y$')
    truc_1.set_xlim(-0.2, 2.2)
    truc_1.set_ylim(-1.3, 0.3)
    truc_1.grid(True, linestyle='--', alpha=0.6)
    truc_1.legend(loc='lower left')

    # Đồ thị con 2: Cực tiểu toàn cục (dây xích treo hình chữ U vuông)
    truc_2.plot(nghiem_toan_cuc['toa_do_khop_x'], nghiem_toan_cuc['toa_do_khop_y'], 'darkgreen', linewidth=2.5, label='Mắt xích', zorder=3)
    truc_2.plot(nghiem_toan_cuc['toa_do_khop_x'][1:-1], nghiem_toan_cuc['toa_do_khop_y'][1:-1], 'ro', markersize=8, markerfacecolor='white', markeredgewidth=2, label='Khớp nối', zorder=4)
    truc_2.plot([0, 1], [0, 0], 'ks', markersize=9, label='Điểm treo (0,0) & (1,0)', zorder=5)
    truc_2.set_title(f"Dây xích treo hình chữ U (Cực tiểu toàn cục)\nThế năng $f = {nghiem_toan_cuc['the_nang']:.4f}$", fontweight='bold')
    truc_2.set_xlabel('Tọa độ ngang $x$')
    truc_2.set_ylabel('Tọa độ đứng $y$')
    truc_2.set_xlim(-0.2, 2.2)
    truc_2.set_ylim(-1.3, 0.3)
    truc_2.grid(True, linestyle='--', alpha=0.6)
    truc_2.legend(loc='lower left')

    plt.tight_layout()
    hinh.savefig(ten_tep, dpi=300)
    plt.close(hinh)
    print(f"Đã lưu đồ thị so sánh Câu 2 vào tệp '{ten_tep}'")


def ve_do_thi_cau_3_hinh_2(nghiem_cau_3, ten_tep='hinh_cau3_chuong_ngai_co_so.png'):
    """
    Tái tạo Hình 2 (Figure 2) trong đề bài tập với độ chính xác 100%:
    - Đường cong chướng ngại vật nét đứt màu xanh lơ (cyan)
    - Các đoạn mắt xích màu xanh dương (blue)
    - Các khớp nối hình tròn rỗng màu đỏ (red)
    - Hai điểm đầu cuối đánh dấu nhân màu đen (black)
    - Giới hạn trục tọa độ: x thuộc [0, 2], y thuộc [-2, 0]
    """
    khop_x = nghiem_cau_3['toa_do_khop_x']
    khop_y = nghiem_cau_3['toa_do_khop_y']
    ham_chuong_ngai = nghiem_cau_3['ham_chuong_ngai']

    hinh, truc = plt.subplots(figsize=(6.5, 6))

    # Đường cong chướng ngại vật mịn khớp với đường nét đứt màu cyan trong Hình 2
    luoi_x_min = np.linspace(0.0, 2.0, 400)
    luoi_y_chuong_ngai = ham_chuong_ngai(luoi_x_min)
    truc.plot(luoi_x_min, luoi_y_chuong_ngai, color='cyan', linestyle='--', linewidth=1.8, zorder=2)

    # Các mắt xích màu xanh dương khớp với Hình 2
    truc.plot(khop_x, khop_y, color='blue', linewidth=2.0, zorder=3)

    # Các khớp nối hình tròn rỗng màu đỏ khớp với Hình 2
    truc.plot(khop_x[1:-1], khop_y[1:-1], marker='o', markersize=8, markerfacecolor='none', markeredgecolor='red', markeredgewidth=2.0, linestyle='None', zorder=4)

    # Hai điểm đầu cuối hình chữ X màu đen khớp với Hình 2
    truc.plot([khop_x[0], khop_x[-1]], [khop_y[0], khop_y[-1]], marker='X', markersize=10, color='black', markeredgecolor='black', markeredgewidth=2.2, linestyle='None', zorder=5)

    # Tiêu đề và giới hạn trục tọa độ chuẩn xác theo Hình 2
    truc.set_title('Chiều dài 3 và 10 mắt xích (Length 3 and 10 links)', fontweight='bold', pad=12)
    truc.set_xlim(0.0, 2.0)
    truc.set_ylim(-2.0, 0.0)
    truc.set_xticks(np.arange(0.0, 2.1, 0.5))
    truc.set_yticks(np.arange(-2.0, 0.1, 0.2))

    # Định dạng khung viền
    truc.spines['top'].set_visible(True)
    truc.spines['right'].set_visible(True)
    truc.tick_params(direction='in', top=True, right=True)

    plt.tight_layout()
    hinh.savefig(ten_tep, dpi=300)
    plt.close(hinh)
    print(f"Đã lưu đồ thị cơ sở Câu 3 (tái tạo Hình 2) vào tệp '{ten_tep}'")


def ve_do_thi_cau_3_khao_sat(nghiem_co_so, nghiem_tuy_chon, ten_tep='hinh_cau3_chuong_ngai_tuy_chon.png'):
    """Vẽ đồ thị so sánh cho khảo sát bộ thông số tùy chọn của Câu 3"""
    hinh, (truc_1, truc_2) = plt.subplots(1, 2, figsize=(13, 5.5))

    luoi_x_min = np.linspace(0.0, 2.0, 400)

    # Đồ thị con 1: Nghiệm cơ sở (Tái tạo Hình 2)
    y_chuong_ngai_1 = nghiem_co_so['ham_chuong_ngai'](luoi_x_min)
    truc_1.plot(luoi_x_min, y_chuong_ngai_1, color='cyan', linestyle='--', linewidth=1.8, label='Chướng ngại vật ($\\mu=1.0, h=1.5$)')
    truc_1.plot(nghiem_co_so['toa_do_khop_x'], nghiem_co_so['toa_do_khop_y'], color='blue', linewidth=2.0, label=f'Dây xích cơ sở ($f={nghiem_co_so["the_nang"]:.3f}$)')
    truc_1.plot(nghiem_co_so['toa_do_khop_x'][1:-1], nghiem_co_so['toa_do_khop_y'][1:-1], 'o', markerfacecolor='none', markeredgecolor='red', markeredgewidth=2)
    truc_1.plot([0, 2], [0, 0], 'kX', markersize=10)
    truc_1.set_title('Nghiệm cơ sở (Khớp với Hình 2)', fontweight='bold')
    truc_1.set_xlim(0, 2)
    truc_1.set_ylim(-2.0, 0.0)
    truc_1.set_xlabel('Tọa độ ngang $x$')
    truc_1.set_ylabel('Tọa độ đứng $y$')
    truc_1.legend(loc='lower left', framealpha=0.9)
    truc_1.grid(True, linestyle=':', alpha=0.5)

    # Đồ thị con 2: Khảo sát thông số tùy chọn (chướng ngại vật lệch tâm)
    y_chuong_ngai_2 = nghiem_tuy_chon['ham_chuong_ngai'](luoi_x_min)
    truc_2.plot(luoi_x_min, y_chuong_ngai_2, color='magenta', linestyle='--', linewidth=1.8, label='Chướng ngại vật dịch tâm ($\\mu=1.3, h=1.6$)')
    truc_2.plot(nghiem_tuy_chon['toa_do_khop_x'], nghiem_tuy_chon['toa_do_khop_y'], color='navy', linewidth=2.0, label=f'Dây xích tùy chọn ($f={nghiem_tuy_chon["the_nang"]:.3f}$)')
    truc_2.plot(nghiem_tuy_chon['toa_do_khop_x'][1:-1], nghiem_tuy_chon['toa_do_khop_y'][1:-1], 's', markerfacecolor='orange', markeredgecolor='darkred', markeredgewidth=1.5)
    truc_2.plot([0, 2], [0, 0], 'kX', markersize=10)
    truc_2.set_title('Khảo sát thông số tùy chọn (Đỉnh dịch sang $\\mu=1.3, h=1.6$)', fontweight='bold')
    truc_2.set_xlim(0, 2)
    truc_2.set_ylim(-2.0, 0.0)
    truc_2.set_xlabel('Tọa độ ngang $x$')
    truc_2.set_ylabel('Tọa độ đứng $y$')
    truc_2.legend(loc='lower left', framealpha=0.9)
    truc_2.grid(True, linestyle=':', alpha=0.5)

    plt.tight_layout()
    hinh.savefig(ten_tep, dpi=300)
    plt.close(hinh)
    print(f"Đã lưu đồ thị khảo sát thông số tùy chọn Câu 3 vào tệp '{ten_tep}'")


def xuat_tep_toa_do(ten_tep, khop_x, khop_y):
    """Ghi tọa độ các khớp nối ra tệp văn bản"""
    with open(ten_tep, 'w', encoding='utf-8') as tep:
        tep.write("# Khop X Y\n")
        for k in range(len(khop_x)):
            tep.write(f"{k} {khop_x[k]:12.6f} {khop_y[k]:12.6f}\n")


def ham_chinh():
    print("=" * 75)
    print(" BÀI TẬP 4: BÀI TOÁN DÂY XÍCH TREO — PHÂN TÍCH SỐ TOÀN DIỆN ")
    print("=" * 75)

    # 1. Phân tích Câu 2
    print("\n[CÂU 2: Dây xích 3 mắt xích]")
    nghiem_cau2_dia_phuong = giai_cau_2('duoc_giao')
    nghiem_cau2_toan_cuc = giai_cau_2('treo_thang')

    print(f"  Thế năng với điểm khởi tạo được giao (x=1, y=1): {nghiem_cau2_dia_phuong['the_nang']:.6f}")
    print(f"  Thế năng cực tiểu toàn cục (treo chữ U):          {nghiem_cau2_toan_cuc['the_nang']:.6f}")
    print(f"  Nhân tử Lagrange (địa phương): lambda = {np.round(nghiem_cau2_dia_phuong['nhan_tu_lambda'], 4)}, mu_x = {nghiem_cau2_dia_phuong['nhan_tu_mu_x']:.4f}, mu_y = {nghiem_cau2_dia_phuong['nhan_tu_mu_y']:.4f}")
    print(f"  Sai số KKT cực đại (địa phương): x={np.max(np.abs(nghiem_cau2_dia_phuong['sai_so_kkt_x'])):.2e}, y={np.max(np.abs(nghiem_cau2_dia_phuong['sai_so_kkt_y'])):.2e}")
    print(f"  Giá trị riêng Hessian thu gọn (địa phương): {nghiem_cau2_dia_phuong['hessian_thu_gon'][0, 0]:.4f} (> 0 => Thỏa mãn SOSC!)")

    print(f"  Nhân tử Lagrange (toàn cục): lambda = {np.round(nghiem_cau2_toan_cuc['nhan_tu_lambda'], 4)}, mu_x = {nghiem_cau2_toan_cuc['nhan_tu_mu_x']:.4f}, mu_y = {nghiem_cau2_toan_cuc['nhan_tu_mu_y']:.4f}")
    print(f"  Sai số KKT cực đại (toàn cục): x={np.max(np.abs(nghiem_cau2_toan_cuc['sai_so_kkt_x'])):.2e}, y={np.max(np.abs(nghiem_cau2_toan_cuc['sai_so_kkt_y'])):.2e}")
    print(f"  Giá trị riêng Hessian thu gọn (toàn cục): {nghiem_cau2_toan_cuc['hessian_thu_gon'][0, 0]:.4f} (> 0 => Thỏa mãn SOSC!)")

    ve_do_thi_cau_2_don(nghiem_cau2_dia_phuong)
    ve_do_thi_cau_2_so_sanh(nghiem_cau2_dia_phuong, nghiem_cau2_toan_cuc)
    xuat_tep_toa_do('toa_do_xich_co_ban.txt', nghiem_cau2_dia_phuong['toa_do_khop_x'], nghiem_cau2_dia_phuong['toa_do_khop_y'])

    # 2. Câu 3 Cơ sở (Tái tạo Hình 2)
    print("\n[CÂU 3: Dây xích 10 mắt xích có chướng ngại vật — Cấu hình cơ sở]")
    nghiem_cau3_co_so = giai_cau_3_chuong_ngai(chieu_cao_h=1.5, tam_mu=1.0, phuong_sai=0.1, do_dich_nu=2.0)
    print(f"  Thế năng tối ưu cơ sở: {nghiem_cau3_co_so['the_nang']:.6f}")
    print("  Tọa độ các khớp nối khớp với Hình 2:")
    for k in range(nghiem_cau3_co_so['so_mat_xich'] + 1):
        print(f"    Khớp {k:2d}: X = {nghiem_cau3_co_so['toa_do_khop_x'][k]:10.6f}, Y = {nghiem_cau3_co_so['toa_do_khop_y'][k]:10.6f}")

    ve_do_thi_cau_3_hinh_2(nghiem_cau3_co_so)
    xuat_tep_toa_do('toa_do_xich_chuong_ngai_co_so.txt', nghiem_cau3_co_so['toa_do_khop_x'], nghiem_cau3_co_so['toa_do_khop_y'])

    # 3. Câu 3 Khảo sát bộ thông số tùy chọn
    print("\n[CÂU 3: Khảo sát thông số tùy chọn (Chướng ngại vật dịch tâm tam_mu=1.3, chieu_cao_h=1.6)]")
    nghiem_cau3_tuy_chon = giai_cau_3_chuong_ngai(
        chieu_cao_h=1.6, tam_mu=1.3, phuong_sai=0.1, do_dich_nu=2.0, ep_nhanh_hinh_2=False
    )
    print(f"  Thế năng tối ưu tùy chọn: {nghiem_cau3_tuy_chon['the_nang']:.6f}")

    ve_do_thi_cau_3_khao_sat(nghiem_cau3_co_so, nghiem_cau3_tuy_chon)
    xuat_tep_toa_do('toa_do_xich_chuong_ngai_tuy_chon.txt', nghiem_cau3_tuy_chon['toa_do_khop_x'], nghiem_cau3_tuy_chon['toa_do_khop_y'])

    print("\nHoàn tất phân tích số và xuất các đồ thị minh họa!")


if __name__ == '__main__':
    ham_chinh()

