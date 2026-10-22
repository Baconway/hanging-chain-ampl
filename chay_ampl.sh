#!/usr/bin/env bash
# ==============================================================================
# Kịch bản hỗ trợ: chay_ampl.sh
# Thực thi các mô hình AMPL cho Bài tập 4 hoặc chạy bộ giải Python đi kèm
# ==============================================================================

set -e

echo "=================================================================="
echo "      BÀI TẬP 4: KỊCH BẢN THỰC THI BÀI TOÁN DÂY XÍCH TREO         "
echo "=================================================================="

if command -v ampl &> /dev/null; then
    echo "Đã phát hiện trình thực thi AMPL: $(command -v ampl)"
    echo ""
    echo "--- [1/2] Đang chạy Câu 2 (xich_co_ban.run) ---"
    ampl xich_co_ban.run
    echo ""
    echo "--- [2/2] Đang chạy Câu 3 (xich_chuong_ngai_vat.run) ---"
    ampl xich_chuong_ngai_vat.run
    echo ""
    echo "Thực thi AMPL hoàn tất thành công."
else
    echo "Lưu ý: Không tìm thấy tệp thực thi 'ampl' trong biến môi trường PATH."
    echo "Để chạy bằng AMPL, hãy cài đặt AMPL và thêm vào PATH,"
    echo "hoặc mở trực tiếp trong VS Code với tiện ích mở rộng AMPL chính thức."
    echo ""
    echo "Đang chạy bộ giải số và vẽ đồ thị Python đi kèm (mo_phong_day_xich.py)..."
    python3 mo_phong_day_xich.py
fi

