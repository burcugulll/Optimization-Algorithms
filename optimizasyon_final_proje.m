veri_tam = load('seeds_dataset.txt');%veri setini yükleme
veri=veri_tam(:,1:7);%veri setini son sütundaki kümeleme değerleri alınmadan yazıldı
N = 50;%parcacık sayısı
D=21;%boyut sayısı
maxiter = 1000;%itearsyon sayısı
CR = 0.9;%dea alg kendine özgü param(çaprazlama oranı)her bir bireyin kromozomozomunu güncellemek için kullanılan parametredir. CR=1, her iki birey arasında tam bir çaprazlama gerçekleşir.CR=0, çaprazlama işlemi gerçekleşmez ve yeni birey tamamen eski birey olarak kalır.
F = 0.5;%dea alg kendine özgü param(mutasyon faktörü)Mutasyon faktörü, popülasyon içindeki bireyler arasında gerçekleştirilen mutasyonun büyüklüğünü kontrol eden bir parametredir F değeri, iki rastgele seçilen birey arasındaki fark vektörüne çarpılarak mutasyon büyüklüğünü belirler. Büyük F değerleri, genetik çeşitliliği artırabilir, ancak aşırı büyük değerler algoritmanın konverjansını etkileyebilir.
alt=0; 
ust=1;
ub=ones(1,D)*ust; 
lb=ones(1,D)*alt;

% Veriyi normalize etme
[rows, cols]=size(veri);%veri setinde satır sutun sayısı döndürme
normalize_veri=zeros(rows,cols);%adlı bir matris oluşturulur ve bu matris, normalize edilmiş veriyi içerecek şekilde boyutlandırılır. Boyutları, orijinal veri matrisi olan veri matrisi ile aynıdır.
for col=1:cols%her sutun normalize edilmesi 
     min_deger=min(veri(:,col));%her sütunun min max değerleri bulundu
    max_deger=max(veri(:,col));
    normalize_veri(:,col)=(veri(:,col)-min_deger)./(max_deger-min_deger);%orijinal değerlerden minimum değeri çıkarıp, minimum ve maksimum arasındaki farka bölerek elde edilir. Normalize edilmiş değer [0, 1] aralığına ölçeklendirilildi
end

% Başlangıç popülasyonunu arama uzayı(0-1) sınırları içinde oluştur
G=rand(N,D);%alt 0 ust 1 old için yazılmadı

for i=1:N
 [ObjVal(i),tahminEdilenKumeMerkezi(:,i)] = objfunc_final(G(i,:), normalize_veri);%pop her bireyin uygunluk değeri ve küme merkezleri hesaplanır
end

iter = 1;
while(iter <= maxiter)
    for i = 1:N
        r = randi(N, 1, 3);%rastgele üç birey seç
        while((i == r(1)) || (i == r(2)) || (i == r(3)) || (r(1) == r(2)) || (r(1) == r(3)) || (r(2) == r(3)))
            r = randi(N, 1, 3);%rastgele 3 birey seçildi ve bu 3 komşu kromozomlar birbirlerinden ve parçaıktan farklı olmalı
        end
        
        % Mutasyon
        Fark = G(r(2), :) - G(r(3), :);%2. ve 3.birey arasındaki fark alındı
        Fark = Fark * F;%mutasyonla çarpıldı
        ind = find(Fark < lb);%sınır aşımı kontrolu
        Fark(ind) = lb(ind);
        ind = find(Fark > ub);
        Fark(ind) = ub(ind);
        
        Gyeni = G(r(1), :) + Fark;%1.birey eklenerek yeni birey oluşturuldu
        ind = find(Gyeni < lb);
        Gyeni(ind) = lb(ind);
        ind = find(Gyeni > ub);
        Gyeni(ind) = ub(ind);
        
        % Çaprazlama
        NewKro = zeros(1, D);%çaprazlama sonucundaki yeni birey için yer açıldı
        sayac = 0;%kaç adet gen çaprazlama ile değşti onu sayıyor
        for j = 1:D
            a = rand();
            if(a < CR)
                NewKro(j) = Gyeni(j);%Döngü, her bir gen için, rastgele bir sayı üretilerek CR değerinden küçükse, o gen çaprazlama ile değiştirilir 
                sayac = sayac + 1;%Yeni bireyin j'inci geni, Gyeni vektöründen alınarak NewKro'ya atanır.
            else
                NewKro(j) = G(i, j);%sayac=0eğer hçbir gen çaprazlama ile değişmediyse rastgele seçilen bir gen, Gyeni vektöründen alınarak değiştirilir. Bu, çeşitliliği artırmak için yapılır.
            end%Yeni bireyin j'inci geni, orijinal popülasyondaki bireyin aynı sıradaki geninden alınarak NewKro'ya atanır.
        end
        
        % Hiçbir gen yeni oluşandan alınmamışsa rassal bir boyut yeniden gelecek şekilde değiştirildi
        if(sayac == 0)%0.1 ihtimalle gyeniden faydalanılmadıEğer hiçbir gen çaprazlama ile değiştirilmediyse
            a = fix(rand() * D) + 1;% Rastgele bir gen seçilir.
            NewKro(a) = Gyeni(a);%Yeni bireyin seçilen rastgele geni, Gyeni vektöründen alınarak değiştirilir.
        end
        
        % Çaprazlama
       % Fitness = feval(objfunc_final, NewKro);
        [Fitness,kumeMerkezi] = objfunc_final(NewKro, normalize_veri);

        if(Fitness < ObjVal(i))% yeni bireyin uygunluk değerinin, mevcut bireyin uygunluk değerinden daha iyi olup olmadığını kontrol eder. eğer yeni birey daha iyi bir çözüm sunuyorsa, bu blok çalıştırılır.
            G(i, :) = NewKro;%mevcut bireyin genetik matrisini (G), yeni oluşturulan bireyin genetik matrisi (NewKro) ile günceller. 
            ObjVal(i) = Fitness;%mevcut bireyin uygunluk değerini (ObjVal(i)) yeni bireyin uygunluk değeri ile günceller.
            tahminEdilenKumeMerkezi(:,i)=kumeMerkezi;%mevcut bireyin tahmin ettiği küme merkezlerini (tahminEdilenKumeMerkezi) günceller. 
        end
    end
    
    GlobalMin = min(ObjVal);%ObjVal vektöründeki uygunluk değerlerinden en küçük olanını bulur.
    fprintf('İter=%d ObjVal=%g\n', iter, GlobalMin);
    if(GlobalMin == 0)
        break;
    end
        iter = iter + 1;

end
