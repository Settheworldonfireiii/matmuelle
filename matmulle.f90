program matmulle
implicit none
real, dimension(4,4) :: A1, A2, A3, A4, A5, A6, R12, R123, R1_4, R1_5, R1_6,TB, RES
real(4)    :: theta1, alpha1, d1, a_1, theta2, alpha2, d2, a_2, &
    theta3, alpha3, d3, a_3, theta4, alpha4, d4, a_4, theta5, alpha5, d5, a_5, &
    theta6, alpha6, d6, a_6
integer :: i

    
    READ(*,*) theta1, alpha1, d1, a_1, theta2, alpha2, d2, a_2, &
        theta3, alpha3, d3, a_3, theta4, alpha4, d4, a_4, &
        theta5, alpha5, d5, a_5, theta6, alpha6, d6, a_6
    A1 = reshape([cos(theta1), sin(theta1), 0., 0., -sin(theta1)*cos(alpha1), &
    cos(theta1)*cos(alpha1), sin(alpha1), 0. , sin(theta1)*sin(alpha1), &
    -cos(theta1)*sin(alpha1), cos(alpha1),0., &
        a_1*cos(theta1), a_1*sin(theta1), d1, 1.], [4, 4])
    A2 = reshape([cos(theta2), sin(theta2), 0., 0., -sin(theta2)*cos(alpha2), &
    cos(theta2)*cos(alpha2), sin(alpha2), 0. , sin(theta2)*sin(alpha2), &
    -cos(theta2)*sin(alpha2), cos(alpha2),0., &
        a_2*cos(theta2), a_2*sin(theta2), d2, 1.], [4, 4])
    A3 = reshape([cos(theta3), sin(theta3), 0., 0., -sin(theta3)*cos(alpha3), &
    cos(theta3)*cos(alpha3), sin(alpha3), 0. , sin(theta3)*sin(alpha3), &
    -cos(theta3)*sin(alpha3), cos(alpha3),0., &
        a_3*cos(theta3), a_3*sin(theta3), d3, 1.], [4, 4])
    A4 = reshape([cos(theta4), sin(theta4), 0., 0., -sin(theta4)*cos(alpha4), &
    cos(theta4)*cos(alpha4), sin(alpha4), 0. , sin(theta4)*sin(alpha4), &
    -cos(theta4)*sin(alpha4), cos(alpha4),0., &
        a_4*cos(theta4), a_4*sin(theta4), d4, 1.], [4, 4])
    A5 = reshape([cos(theta5), sin(theta5), 0., 0., -sin(theta5)*cos(alpha5), &
    cos(theta5)*cos(alpha5), sin(alpha5), 0. , sin(theta5)*sin(alpha5), &
    -cos(theta5)*sin(alpha5), cos(alpha5),0., &
        a_5*cos(theta5), a_5*sin(theta5), d5, 1.], [4, 4])
    A6 = reshape([cos(theta6), sin(theta6), 0., 0., -sin(theta6)*cos(alpha6), &
    cos(theta6)*cos(alpha6), sin(alpha6), 0. , sin(theta6)*sin(alpha6), &
    -cos(theta6)*sin(alpha6), cos(alpha6),0., &
        a_6*cos(theta6), a_6*sin(theta6), d6, 1.], [4, 4])
    TB = reshape([-1.,0.,0.,0.,0.,-1.,0.,0.,0.,0.,1.,0.,-5.,0.,-5.,1.], [4, 4])
    R12 = matmul(A1, A2)
    R123 = matmul(R12, A3)
    R1_4 = matmul(R123, A4)
    R1_5 = matmul(R1_4, A5)
    R1_6 = matmul(R1_5, A6)
    RES = matmul(TB, R1_4)
    do i=1,4
        Print *, A1(i,:)
    end do
    do i=1,4
        Print *, A2(i,:)
    end do
    do i=1,4
        Print *, A3(i,:)
    end do
    do i=1,4
        Print *, A4(i,:)
    end do
    do i=1,4
        Print *, A5(i,:)
    end do
    do i=1,4
        Print *, A6(i,:)
    end do
    do i=1,4
        Print *, R12(i,:)
    end do
    do i=1,4
        Print *, R123(i,:)
    end do
    do i=1,4
        Print *, R1_4(i,:)
    end do
    do i=1,4
        Print *, R1_5(i,:)
    end do
    do i=1,4
        Print *, R1_6(i,:)
    end do
    do i=1,4
        Print *, RES(i,:)
    end do
end program matmulle
