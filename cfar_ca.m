function [cfar_out,cfar_ref] = cfar_ca(input_array,pfa,ref_num,pro_num)
    % pro_num ref_num都是单边单元个数
    % 为设置边缘保护，存在边缘漏检可能性
    % input_array为复矩阵

    alpha = 2*ref_num*(pfa^(-1/2/ref_num)-1);
    % alpha = 14;
    [row,col] = size(input_array);
    cfar_in = abs(input_array).^2;
    cfar_out = zeros(row,col);
    cfar_ref = cfar_out;
    for i = 1:row
        for j = ref_num+pro_num+1:col-(ref_num+pro_num)
            cfar_ref(i,j) = 1/2/ref_num*alpha*(sum(cfar_in(i,j-(ref_num+pro_num):j-pro_num-1))+sum(cfar_in(i,j+pro_num+1:j+pro_num+ref_num)));
            cfar_out(i,j) = cfar_in(i,j)>cfar_ref(i,j);
        end
    end
end